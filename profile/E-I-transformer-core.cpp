// -*- coding: utf-8 -*-
/* 
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and
 * library for computing electromagnetic fields.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dirichlet_condition.hpp"
#include "grid_compiler.hpp"
#include "grid_patch.hpp"
#include "linear_vertex_distribution.hpp"
#include "material.hpp"
#include "newton.hpp"
#include "reluctivity_models.hpp"
#include "shape.hpp"

using namespace PACKAGE_NS;
using boost::units::si::ampere;
using boost::units::si::ampere_per_metre;
using boost::units::si::ampere_per_metre_squared;
using boost::units::si::metre;
using boost::units::si::metre_per_henry;
using boost::units::si::radians;
using boost::units::si::second;
using boost::units::si::siemens_per_metre;
using boost::units::si::square_metre;
using boost::units::si::tesla;
using boost::units::si::weber_per_metre;
using boost::units::si::constants::codata::mu_0;

constexpr double kPi = boost::math::double_constants::pi;

// Geometry
auto air_surround = 2e-2*metre;
auto coil_dim_x = 1e-2*metre;
auto coil_dim_y = 4e-2*metre;
auto iron_width = 2e-2*metre;
auto coil_space_top = 5e-3*metre;
auto coil_space_left_side = 3e-3*metre;
auto coil_space_right_side = 7e-3*metre;
auto domain_height = coil_dim_y+coil_space_top+iron_width+air_surround;
auto domain_width = (2.*iron_width+coil_dim_x+coil_space_left_side
		     +coil_space_right_side+air_surround);
auto grid_division = 18e2/metre;

// Electro - magnetic
auto current_density = 2.5e6*ampere_per_metre_squared;
auto frequency = 50./second;

struct M30035AH {
   static constexpr std::size_t kNumberOfValues = 19;
   static const MagneticFieldIntensityReal kFieldIntensity[kNumberOfValues];
   static const MagneticPolarizationReal kPolarization[kNumberOfValues];
};

const MagneticFieldIntensityReal M30035AH::kFieldIntensity [] = {
   .0 * ampere_per_metre, 31.4 * ampere_per_metre,
   41.4 * ampere_per_metre, 48.2 * ampere_per_metre,
   54.3 * ampere_per_metre, 60.4 * ampere_per_metre,
   67.1 * ampere_per_metre, 74.9 * ampere_per_metre,
   84.2 * ampere_per_metre, 96.3 * ampere_per_metre,
   113. * ampere_per_metre, 137. * ampere_per_metre,
   179. * ampere_per_metre, 266. * ampere_per_metre,
   521. * ampere_per_metre, 1380. * ampere_per_metre,
   3400. * ampere_per_metre, 6610. * ampere_per_metre,
   11100. * ampere_per_metre
};

const MagneticPolarizationReal M30035AH::kPolarization [] = {
   .0 * tesla, .1 * tesla, .2 * tesla, .3 * tesla, .4 * tesla,
   .5 * tesla, .6 * tesla, .7 * tesla, .8 * tesla, .9 * tesla,
   1. * tesla, 1.1 * tesla, 1.2 * tesla, 1.3 * tesla,
   1.4 * tesla, 1.5 * tesla, 1.6 * tesla, 1.7 * tesla, 1.8 * tesla
};

auto iron_fill_factor = .97;

auto material = Material::New(
   SplineReluctivity::New(
      M30035AH::kFieldIntensity,
      M30035AH::kFieldIntensity+M30035AH::kNumberOfValues,
      M30035AH::kPolarization));


int main()
{
   using CPoint = CartesianSystem::Point;
   using CRectangle = Rectangle<CartesianSystem>;
   using GridPatchCartesian =
      GridPatch<CartesianSystem, AnalysisType >;
   using CurrentDensity = AnalysisType::CurrentDensity;
   using Time = AnalysisType::Time;
   
   auto vd_x = LinearVertexDistribution::New(domain_width*grid_division+1);
   auto vd_y = LinearVertexDistribution::New(domain_height*grid_division+1);

   auto p0 = CPoint(.0*metre, .0*metre);
   auto p1 = CPoint(domain_width, domain_height);

   auto grid = GridPatchCartesian::New(p0, p1, vd_x, vd_y);

   auto coil_llc = CPoint(iron_width+coil_space_left_side, .0*metre);
   grid->SetLoad(
      CRectangle::New(coil_llc, coil_dim_x, coil_dim_y),
      [](Time t) -> CurrentDensity {
	 return current_density*sin(2*kPi*(frequency*t));
      });

   auto iron_height_in = coil_dim_y + coil_space_top;
   auto iron_width_total = (2.*iron_width+ coil_space_left_side
			    +coil_space_right_side+coil_dim_x);
   ConstShapePointer<CartesianSystem> iron_rects[] = {
      CRectangle::New(
	 CPoint(.0*metre, .0*metre), iron_width, iron_height_in),
      CRectangle::New(
	 CPoint(
	    iron_width+coil_dim_x+coil_space_left_side
	    +coil_space_right_side, .0*metre), iron_width, iron_height_in),
      CRectangle::New(
	 CPoint(.0*metre, iron_height_in), iron_width_total, iron_width)
   };

   for(std::size_t i=2; i; --i)
      grid->SetMaterial(iron_rects[i], material);

   auto dbc = DirichletCondition<AnalysisType>::New();

   grid->SetBoundaryConditionWest(dbc);
   grid->SetBoundaryConditionNorth(dbc);
   grid->SetBoundaryConditionEast(dbc);

   fi::GridCompiler<AnalysisType> grid_compiler(grid);

   const std::size_t number_of_timesteps = 16;
   AnalysisType::Time time_step = 1./number_of_timesteps/frequency;

   auto grid_graph = grid_compiler(time_step);
   
   for(std::size_t k=0; k < number_of_timesteps; ++k) {
      std::size_t max_iter_steps = 24;
      Newton(grid_graph, 1e-5, max_iter_steps);
      NextTimeStep(grid_graph);
   }
}

/*
 * Local Variables: 
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
