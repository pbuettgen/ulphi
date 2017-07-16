// -*- coding: utf-8 -*-
/*
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and library for computing
 * electromagnetic fields.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*!
 *  \file
 *
 *  \brief (Inhomogenious) neumann boundary condition
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef NEUMANN_VERTEX_HPP_vsz6TMPa_
#define NEUMANN_VERTEX_HPP_vsz6TMPa_

#include "config.h"

#include "grid_graph-fwd.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Vertex type implementing the Neumann boundary condition
      /*!
       *  The Neumann boundary condition forces the first derivation to a
       *  fixed value.  In the case considered here, the first derivation is
       *  nothing else then the magnetic flux density.
       *
       *  \tparam Analysis defines properties related to the selected
       *  analysis, especially the data type for magnetic quantities.
       *
       *  \tparam Traits defines the in detail behaviour of this class.
       */
      template <typename Analysis, template <typename> class Traits>
      struct NeumannVertex:
         StandardVertex<Analysis, Traits>,
         Traits<Analysis>::NeumannVertexProperties
      {
      protected:
         //! Base class type
         using Base = StandardVertex<Analysis, Traits>;
         using VertexTraits = Traits<Analysis>;
         using NeumannVertexProperties =
            typename VertexTraits::NeumannVertexProperties;
         using GridGraphPtr = typename Base::GridGraphPtr;
         using AssemblingMatrix = typename Base::AssemblingMatrix;
         using CoeffCoeffDerivedPair = typename Base::CoeffCoeffDerivedPair;
         using CurrentDensity = typename Analysis::CurrentDensity;
         using MagneticVectorPotential = typename Base::MagneticVectorPotential;
         using Point = typename Base::Point;
         using RhsVector = typename Base::RhsVector;
         using VertexDescriptor = typename Base::VertexDescriptor;
         using LoadDefPointer = ConstLoadDefinitionPointer<Analysis>;
         //! Magnetic flux density type
         using MagneticFluxDensity = typename Analysis::MagneticFluxDensity;
         using ReluctivityPassType = PassType<MagneticReluctivity>;
         using FluxDensityPassType = PassType<MagneticFluxDensity>;
         using LengthPassType = PassType<Length>;
         using EdgeProps = typename Base::EdgeProps;

      public:
         //! Default initialization
         NeumannVertex() : Base(), NeumannVertexProperties()
         {
         }

         //! Initialize a Neumann vertex
         /*!
          * @param position is the vertex' position within the global
          * cartesian coordinate system
          *
          * @param reluctivity_outside is the homogeneous reluctivity
          *  outside the computation domain
          *
          * @param flux_density_outside is the magnetic flux density
          *  outside the computation domain
          *
          * @param delta is the grid's edge length perpendicular to
          * the boundary
          *
          * @param parallelity_factor defines weather the boundary is
          * parallel or antiparallel to the coordinate system's axis.
          *
          * @param args are forwarded to the base class
          */
         template <typename ... Args>
         NeumannVertex(
            const Point& position,
            MagneticReluctivity reluctivity_outside,
            PassType<MagneticFluxDensity> flux_density_outside, Length delta,
            Parallelity parallelity_factor, Args ... args)
            : Base(position, args ...),
              NeumannVertexProperties(
                 reluctivity_outside, flux_density_outside, delta,
                 parallelity_factor)
         {
         }

         //! Destruction
         virtual ~NeumannVertex()
         {
         }

         void AddToEquationSystem(
            VertexDescriptor vertex_id, GridGraphPtr grid_graph,
            AssemblingMatrix& system_triplets, RhsVector& rhs) const override
         {
            rhs(this->equation_number()) += in_base_units(this->NeumannRhs());

            Base::AddToEquationSystem(
               vertex_id, grid_graph, system_triplets, rhs);
         }

         void AddToEquationSystem(
            VertexDescriptor vertex_id, GridGraphPtr grid_graph,
            AssemblingMatrix& system_triplets,
            AssemblingMatrix& jacobi_triplets, RhsVector& rhs) const override
         {
            rhs(this->equation_number()) += in_base_units(this->NeumannRhs());

            Base::AddToEquationSystem(
               vertex_id, grid_graph, system_triplets, jacobi_triplets, rhs);
         }
      };
   }
}

#endif // NEUMANN_VERTEX_HPP_vsz6TMPa_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
