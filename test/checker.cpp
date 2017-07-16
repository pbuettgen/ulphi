// -*- coding: utf-8 -*-
/*
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and library for
 * computing electromagnetic fields.
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

#include <cmath>
#include <limits>
#include <iostream>
#include <tuple>
#include <random>

#include "config.h"

#include <boost/units/cmath.hpp>
#include <boost/units/io.hpp>
#include <boost/units/systems/si.hpp>

#include <gtest/gtest.h>

#include "analysis.hpp"
#include "circle.hpp"
#include "coordinate_systems.hpp"
#include "dirichlet_condition.hpp"
#include "exp_vertex_distribution.hpp"
#include "grid_assembly.hpp"
#include "grid_compiler.hpp"
#include "grid_graph.hpp"
#include "grid_patch.hpp"
#include "linear_vertex_distribution.hpp"
#include "pow.hpp"
#include "material.hpp"
#include "reluctivity_models.hpp"
#include "rectangle.hpp"
#include "standard_vertex.hpp"
#include "types.hpp"

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


struct MagCurvChecks
{
   template <typename HIterator, typename BiIterator>
   static double CompError(
      const ScalarReluctivity& reluctivity,
      HIterator h_begin, HIterator h_end, BiIterator bi_begin)
   {
      double error = .0;

      for (; h_begin != h_end; ++h_begin, ++bi_begin)
      {
         if ((.0*ampere_per_metre == *h_begin) or (.0*tesla == *bi_begin))
            continue;
         const auto b = mu_0 * (*h_begin) + (*bi_begin);
         const auto nu_should_be = (*h_begin) / b;
         const auto reluctivity_reluctivity_derived = reluctivity(b);
         const auto nu = get<0>(reluctivity_reluctivity_derived);
         error += Pow<2>(in_base_units(nu - nu_should_be));
      }

      return sqrt(error);
   }

   static void CheckDeriv(const ScalarReluctivity& reluctivity)
   {
      for (MagneticFluxDensityReal b = .075*tesla; b < 1.625*tesla; b+=.075*tesla)
      {
         const MagneticFluxDensityReal delta_b
            = 2*sqrt(std::numeric_limits<
                     MagneticFluxDensityReal::value_type>::epsilon())*b;
         MagneticFluxDensityReal b_plus_delta_b = b + delta_b;
         const MagneticFluxDensityReal db = b_plus_delta_b - b;

         const auto rrd1 = reluctivity(b-delta_b/2.);
         const auto rrd2 = reluctivity(b+delta_b/2.);
         const auto rrd = reluctivity(b);

         const double secant
            = in_base_units((get<0>(rrd2)-get<0>(rrd1))/db);
         const double av_deriv
            = in_base_units(get<1>(rrd));

         EXPECT_LE(abs(secant-av_deriv)/av_deriv, 8e-2)
            << "@ " << b << " Secant: " << secant << " Deriv: " << av_deriv;
      }
   }

   static void CheckReluctAtBZero(const ScalarReluctivity& reluctivity)
   {
      const auto reluctivity0 = reluctivity(.0*tesla);
      EXPECT_TRUE(isfinite(get<0>(reluctivity0)));
      EXPECT_TRUE(isfinite(get<1>(reluctivity0)));
   }
};

using AnalysisTypes = testing::Types<TimeTransientAnalysis<>,
                                     MagnetoHarmonicAnalysis<> >;

template <typename Analysis, typename CoordinateSystem>
struct AnalysisAndCSys {
   using AnalysisType = Analysis;
   using CoordinateSystemType = CoordinateSystem;
};

using AnalysisAndCSysTypes = testing::Types<
   AnalysisAndCSys<TimeTransientAnalysis<>, CartesianSystem>,
   AnalysisAndCSys<TimeTransientAnalysis<>, PolarSystem>,
   AnalysisAndCSys<MagnetoHarmonicAnalysis<>, CartesianSystem>,
   AnalysisAndCSys<MagnetoHarmonicAnalysis<>, PolarSystem> >;

// --- ConstantReluctivity ----------------------------------------
struct ConstantReluctivityTest : testing::Test {
   static constexpr double RelTestPermeability()
   {
      return 1e3;
   }

   static Permeability TestPermeability()
   {
      return RelTestPermeability()*mu_0;
   }

   static MagneticReluctivity TestReluctivity()
   {
      return 1./TestPermeability();
   }
};

TEST_F(ConstantReluctivityTest, InitializeByRelPermeability) {
   ConstantReluctivity constant_reluctivity(RelTestPermeability());
   const auto cr_eval = constant_reluctivity(1.*tesla);
   EXPECT_FLOAT_EQ(in_base_units(get<0>(cr_eval)),
                   in_base_units(TestReluctivity()));
   EXPECT_FLOAT_EQ(in_base_units(get<1>(cr_eval)), .0);
}

TEST_F(ConstantReluctivityTest, InitializeByPermeability) {
   ConstantReluctivity constant_reluctivity(TestPermeability());
   const auto cr_eval = constant_reluctivity(1.*tesla);
   EXPECT_FLOAT_EQ(in_base_units(get<0>(cr_eval)),
                   in_base_units(TestReluctivity()));
   EXPECT_FLOAT_EQ(in_base_units(get<1>(cr_eval)), .0);
}

TEST_F(ConstantReluctivityTest, InitilizeByReluctivity) {
   ConstantReluctivity constant_reluctivity(TestReluctivity());
   const auto cr_eval = constant_reluctivity(1.*tesla);
   EXPECT_FLOAT_EQ(in_base_units(get<0>(cr_eval)),
                   in_base_units(TestReluctivity()));
   EXPECT_FLOAT_EQ(in_base_units(get<1>(cr_eval)), .0);
}

// --- ArctanReluctivity ------------------------------------------
struct ArctanReluctivityTest
   : M30035AH, MagCurvChecks, testing::Test {
   ArctanReluctivityTest()
      : arctan_reluctivity_(
         kFieldIntensity, kFieldIntensity+kNumberOfValues, kPolarization)
   {
   }

   // expected values
   static constexpr double kMurShouldBe = 7556;
   static MagneticPolarizationReal JsShouldBe()
   {
      return 1.74 * tesla;
   }
   static constexpr double MaxError() {
      return 5e3;
   }

   ArctanReluctivity arctan_reluctivity_;
};

TEST_F(ArctanReluctivityTest, InitialMuR) {
   EXPECT_FLOAT_EQ(
      round(arctan_reluctivity_.initial_mu_r()), kMurShouldBe);
}

TEST_F(ArctanReluctivityTest, SaturationPolarization) {
   EXPECT_FLOAT_EQ(
      round(100.*in_base_units(
               arctan_reluctivity_.saturation_polarization()))/100.,
      in_base_units(JsShouldBe()));
}

TEST_F(ArctanReluctivityTest, DerivationFromData) {
   const auto error = CompError(
      arctan_reluctivity_,
      kFieldIntensity, kFieldIntensity+kNumberOfValues, kPolarization);

   EXPECT_LE(error, MaxError());
}

TEST_F(ArctanReluctivityTest, Derivative)
{
   CheckDeriv(arctan_reluctivity_);
}

TEST_F(ArctanReluctivityTest, ReluctivityAtBZero)
{
   CheckReluctAtBZero(arctan_reluctivity_);
}

// --- SqrtReluctivity --------------------------------------------
struct SqrtReluctivityTest
   : M30035AH, MagCurvChecks, testing::Test
{
   SqrtReluctivityTest()
      : sqrt_reluctivity_(
         kFieldIntensity, kFieldIntensity+kNumberOfValues, kPolarization)
   {
   }

   static constexpr double kMurShouldBe = 7538;
   static MagneticPolarizationReal JsShouldBe()
   {
      return 1.64 * tesla;
   }
   static constexpr double kAShouldBe = 0.152;
   static constexpr double MaxError() {
      return 5e3;
   }

   SqrtReluctivity sqrt_reluctivity_;
};

TEST_F(SqrtReluctivityTest, InitialMuR)
{
   EXPECT_FLOAT_EQ(round(sqrt_reluctivity_.initial_mu_r()), kMurShouldBe);
}

TEST_F(SqrtReluctivityTest, SaturationPolarization)
{
   EXPECT_FLOAT_EQ(
      round(100.*in_base_units(
               sqrt_reluctivity_.saturation_polarization()))/100.,
      in_base_units(JsShouldBe()));
}

TEST_F(SqrtReluctivityTest, ParameterA)
{
   EXPECT_FLOAT_EQ(round(1e3*sqrt_reluctivity_.a())/1e3, kAShouldBe);
}

TEST_F(SqrtReluctivityTest, DerivationFromData)
{
   const double error = CompError(
      sqrt_reluctivity_,
      kFieldIntensity, kFieldIntensity+kNumberOfValues, kPolarization);

   EXPECT_LE(error, MaxError());
}

TEST_F(SqrtReluctivityTest, Derivative)
{
   CheckDeriv(sqrt_reluctivity_);
}

TEST_F(SqrtReluctivityTest, ReluctivityAtBZero)
{
   CheckReluctAtBZero(sqrt_reluctivity_);
}

// --- SplineReluctivity ------------------------------------------
struct SplineReluctivityTest
   : M30035AH, MagCurvChecks, testing::Test {
   SplineReluctivityTest()
      : spline_reluctivity_(
         kFieldIntensity, kFieldIntensity+kNumberOfValues, kPolarization)
   {
   }

   SplineReluctivity spline_reluctivity_;
};

TEST_F(SplineReluctivityTest, DerivationFromData)
{
   const auto error = CompError(
      spline_reluctivity_,
      kFieldIntensity, kFieldIntensity+kNumberOfValues, kPolarization);

   EXPECT_LE(error, 1.);
}

TEST_F(SplineReluctivityTest, ReluctivityAtBZero)
{
   CheckReluctAtBZero(spline_reluctivity_);
}

// --- CombinedReluctivity ----------------------------------------
struct CombinedReluctivityData {
   static constexpr double Mur1()
   {
      return 1.;
   }

   static constexpr double Mur2()
   {
      return .5;
   }

   static constexpr double Kfe()
   {
      return .5;
   }
};

template <typename> struct CombinedReluctivityTestBase;

template<> struct CombinedReluctivityTestBase<FluxParallel>
   : CombinedReluctivityData
{
   static Permeability MuEff()
   {
      return mu_0 * (Kfe()*Mur1() - Kfe()*Mur2() + Mur2());
   }
};

template<>  struct CombinedReluctivityTestBase<FluxPerpendicular>
   : CombinedReluctivityData
{
   static Permeability MuEff()
   {
      return mu_0 * Mur1()*Mur2()/(Kfe()*Mur2() - Kfe()*Mur1() + Mur1());
   }
};

template <typename FluxOrientation>
struct CombinedReluctivityTest
   : CombinedReluctivityTestBase<FluxOrientation>, testing::Test
{
   CombinedReluctivityTest() : combined_reluctivity_(
      ConstantReluctivity::New(this->Mur1()),
      ConstantReluctivity::New(this->Mur2()), this->Kfe())
   {
   }

   CombinedReluctivity<FluxOrientation> combined_reluctivity_;
};

using FluxOrientations = ::testing::Types<FluxParallel, FluxPerpendicular>;
TYPED_TEST_CASE(CombinedReluctivityTest, FluxOrientations);

TYPED_TEST(CombinedReluctivityTest, Reluctivity)
{
   const auto cr_result = this->combined_reluctivity_(1.*tesla);

   EXPECT_FLOAT_EQ(
      in_base_units(get<0>(cr_result)), in_base_units(1./this->MuEff()));
}

// --- Material ---------------------------------------------------
struct MaterialTest : testing::Test
{
   static Permeability Mu1()
   {
      return 1. * mu_0;
   }
   static Permeability Mu2()
   {
      return .5 * mu_0;
   }
   static ElectricalConductivity Sigma1()
   {
      return 10. * siemens_per_metre;
   }

   MaterialTest() : material_(
      _anisotropic_reluctivity = std::make_tuple(
         ConstantReluctivity::New(Mu1()),
         ConstantReluctivity::New(Mu2())),
      _electrical_conductivity = Sigma1())
   {
   }

   Material material_;
};

TEST_F(MaterialTest, ReluctivityFirstDirection)
{
   using std::get;

   EXPECT_FLOAT_EQ(
      in_base_units(
         get<0> (material_.EvaluateReluctivity(1.*tesla, Axis::kFirst))),
      in_base_units(1./Mu1()));
}

TEST_F(MaterialTest, ReluctivitySecondDirection)
{
   using std::get;

   EXPECT_FLOAT_EQ(
      in_base_units(
         get<0>(material_.EvaluateReluctivity(1.*tesla, Axis::kSecond))),
      in_base_units(1./Mu2()));
}

TEST_F(MaterialTest, ElectricalConductivity)
{
   EXPECT_FLOAT_EQ(
      in_base_units(material_.EvaluateConductivity()),
      in_base_units(Sigma1()));
}

// --- LoadDefinition ---------------------------------------------

template <typename Analysis>
struct LoadDefinitionTest : testing::Test {

   using Time = typename Analysis::Time;
   using LoadDef = LoadDefinition<Analysis>;
   using CurrentDensity = typename Analysis::CurrentDensity;
};

TYPED_TEST_CASE(LoadDefinitionTest, AnalysisTypes);

TYPED_TEST(LoadDefinitionTest, InitByConstCurrentDensity)
{
   using CurrentDensity = typename TypeParam::CurrentDensity;
   using Time = typename TypeParam::Time;
   using LoadDef = LoadDefinition<TypeParam>;

   LoadDef load_definition(CurrentDensity::from_value(1.));

   EXPECT_FLOAT_EQ(
      in_base_units(abs(load_definition(Time::from_value(1.)))), 1.);
}

TYPED_TEST(LoadDefinitionTest, InitByVarCurrentDensity)
{
   using Time = typename TypeParam::Time;
   using CurrentDensity = typename TypeParam::CurrentDensity;
   using LoadDef = LoadDefinition<TypeParam>;

   LoadDef load_definition(
      [](Time t) -> CurrentDensity {
         return CurrentDensity::from_value(1.)*t/Time::from_value(1.);
      });

   for(Time t=Time::from_value(.0); t < Time::from_value(10.);
       t+=Time::from_value(1.))
   {
      EXPECT_FLOAT_EQ(
         in_base_units(abs(load_definition(t))), in_base_units(t));
   }
}

TYPED_TEST(LoadDefinitionTest, InitByConstCurrent)
{
   using Time = typename TypeParam::Time;
   using LoadDef = LoadDefinition<TypeParam>;
   using Current = typename TypeParam::Current;

   LoadDef load_definition(
      Current::from_value(1e6), 1.*square_metre);

   EXPECT_FLOAT_EQ(
      in_base_units(abs(load_definition(Time::from_value(1.)))), 1e6);
}

TYPED_TEST(LoadDefinitionTest, InitByVarCurrent)
{
   using Time = typename TypeParam::Time;
   using Current = typename TypeParam::Current;
   using LoadDef = LoadDefinition<TypeParam>;

   LoadDef load_definition(
      [](Time t) -> Current {
         return Current::from_value(1e6)*t/Time::from_value(1.);
      }, 1.*square_metre);

   for(Time t=Time::from_value(.0); t < Time::from_value(10.);
       t+=Time::from_value(1.))
   {
      EXPECT_FLOAT_EQ(
         in_base_units(abs(load_definition(t))), 1e6*in_base_units(t));
   }
}

// --- Circle -----------------------------------------------------

struct CircleRadius {
   static Length Radius()
   {
      return 1.*metre;
   }
};

template <typename> struct CircleTestBase;

template <> struct CircleTestBase<CartesianSystem> : CircleRadius
{
   using Point = typename CartesianSystem::Point;

   Point Origin() {
      return Point(.0*metre, .0*metre);
   }

   Point MakePointInside()
   {
      const auto limit = .999*in_base_units(Radius());
      std::uniform_real_distribution<> urd_r(-limit, limit);
      const auto radius = urd_r(rand_source_);
      std::uniform_real_distribution<> urd_a(-abs(radius), abs(radius));
      const auto a = urd_a(rand_source_);
      const auto b = sqrt(Pow<2>(radius)-Pow<2>(a));

      return Point(a*metre, b*metre);
   }

   Point MakePointOutside()
   {
      using namespace std;

      constexpr auto RandMax
         = numeric_limits<
            typename uniform_real_distribution<>::result_type>::max();

      uniform_real_distribution<> urd_radius(
         1.0001*in_base_units(Radius()), RandMax);
      const auto radius = urd_radius(rand_source_);
      uniform_real_distribution<> urd_a(-abs(radius), abs(radius));
      const auto a = urd_a(rand_source_);
      const auto b = sqrt(Pow<2>(radius)-Pow<2>(a));

      return Point(a*metre, b*metre);
   }

private:
   std::mt19937 rand_source_;
};

template <> struct CircleTestBase<PolarSystem> : CircleRadius
{
   using Point = typename PolarSystem::Point;

   Point Origin() {
      return Point(.0*metre, .0*radians);
   }

   Point MakePointInside()
   {
      std::uniform_real_distribution<> urd_radius(
         .0, .999*in_base_units(Radius()));
      std::uniform_real_distribution<> urd_phi(.0, 2*kPi);

      return Point(
         urd_radius(rand_source_)*metre, urd_phi(rand_source_)*radians);
   }

   Point MakePointOutside()
   {
      using namespace std;

      constexpr auto RandMax = numeric_limits<
         typename uniform_real_distribution<>::result_type>::max();

      uniform_real_distribution<> urd_radius(
         1.001*in_base_units(Radius()), RandMax);
      uniform_real_distribution<> urd_phi(.0, 2*kPi);

      return Point(
         urd_radius(rand_source_)*metre, urd_phi(rand_source_)*radians);
   }

private:
   std::mt19937 rand_source_;
};

template <typename CoordinateSystem>
struct CircleTest
   : CircleTestBase<CoordinateSystem>, testing::Test
{
   CircleTest()
      : circle_(this->Origin(), this->Radius())
   {
   }

   Circle<CoordinateSystem> circle_;
};

using CoordinateSystems = testing::Types<CartesianSystem, PolarSystem>;
TYPED_TEST_CASE(CircleTest, CoordinateSystems);

TYPED_TEST(CircleTest, PointIsInside)
{
   using Point = typename TypeParam::Point;

   std::function<bool (Point)> test_pred = [this](Point p) -> bool {
      return this->circle_.Within(p);
   };

   for (std::size_t count = 20; count>0; --count)
   {
      // EXPECT_PRED1(test_pred, this->MakePointInside());
      EXPECT_TRUE(this->circle_.Within(this->MakePointInside()));
   }
}

TYPED_TEST(CircleTest, PointIsOutside)
{
   using Point = typename TypeParam::Point;

   std::function<bool (Point)> test_pred = [this](Point p) -> bool {
      return not this->circle_.Within(p);
   };

   for (std::size_t count = 20; count>0; --count)
   {
      // EXPECT_PRED1(test_pred, this->MakePointOutside());
      EXPECT_FALSE(this->circle_.Within(this->MakePointOutside()));
   }
}

// --- Rectangle --------------------------------------------------

template <typename CoordinateSystem> struct RectangleTestBase;

template <> struct RectangleTestBase<CartesianSystem>
{
   using Point = CartesianSystem::Point;

   Point LowerLeftCorner()
   {
      return Point(1.*metre, 1.*metre);
   }

   Length XDim()
   {
      return 2.*metre;
   }

   Length YDim()
   {
      return 1.*metre;
   }

   Point MakePointInside()
   {
      using std::get;
      using std::uniform_real_distribution;

      uniform_real_distribution<> urd_x(0, .999*in_base_units(XDim()));
      uniform_real_distribution<> urd_y(0, .999*in_base_units(YDim()));

      return Point(
         get<0>(LowerLeftCorner()) + urd_x(rand_source_)*metre,
         get<1>(LowerLeftCorner()) + urd_y(rand_source_)*metre);
   }

   Point MakePointOutside()
   {
      using namespace std;

      const auto RandMax = numeric_limits<
         typename uniform_real_distribution<>::result_type>::max();

      uniform_real_distribution<> urd_x(
         1.001*in_base_units(get<0>(LowerLeftCorner())+XDim()), RandMax);
      uniform_real_distribution<> urd_y(
         1.001*in_base_units(get<1>(LowerLeftCorner())+YDim()), RandMax);

      return Point(urd_x(rand_source_)*metre, urd_y(rand_source_)*metre);
   }

private:
   std::mt19937 rand_source_;
};

template <> struct RectangleTestBase<PolarSystem>
{
   using Point = PolarSystem::Point;

   Point LowerLeftCorner() {
      return Point(.5*metre, .0*radians);
   }

   Length XDim()
   {
      return .5*metre;
   }

   Angle YDim()
   {
      return kPi/2*radians;
   }

   Point MakePointInside()
   {
      using std::get;
      using std::uniform_real_distribution;

      uniform_real_distribution<> urd_x(0, .999*in_base_units(XDim()));
      uniform_real_distribution<> urd_y(0, .999*in_base_units(YDim()));

      return Point(
         get<0>(LowerLeftCorner()) + urd_x(rand_source_)*metre,
         get<1>(LowerLeftCorner()) + urd_y(rand_source_)*radians);
   }

   Point MakePointOutside()
   {
      using namespace std;

      const auto RandMax = numeric_limits<
         typename uniform_real_distribution<>::result_type>::max();

      uniform_real_distribution<> urd_x(
         1.001*in_base_units(get<0>(LowerLeftCorner())+XDim()), RandMax);
      uniform_real_distribution<> urd_y(
         1.001*in_base_units(get<1>(LowerLeftCorner())+YDim()), RandMax);

      return Point(urd_x(rand_source_)*metre, urd_y(rand_source_)*radians);
   }

private:
   std::mt19937 rand_source_;
};

template <typename CoordinateSystem>
struct RectangleTest
   : RectangleTestBase<CoordinateSystem>, testing::Test
{
   RectangleTest()
      : rectangle_(this->LowerLeftCorner(), this->XDim(), this->YDim())
   {
   }

   Rectangle<CoordinateSystem> rectangle_;
};

TYPED_TEST_CASE(RectangleTest, CoordinateSystems);

TYPED_TEST(RectangleTest, PointIsInside)
{
   for (unsigned count = 20; count; --count)
   {
      EXPECT_TRUE(this->rectangle_.Within(this->MakePointInside()));
   }
}

TYPED_TEST(RectangleTest, PointIsOutside)
{
   for (unsigned count = 20; count; --count)
   {
      EXPECT_FALSE(this->rectangle_.Within(this->MakePointOutside()));
   }
}

// --- LinearVertexDistribution -----------------------------------
struct LinearVertexDistributionTest : testing::Test
{
   static constexpr std::size_t NumberOfVertices() {
      return 4;
   }

   LinearVertexDistributionTest()
      : linear_vertex_distribution_(NumberOfVertices())
   {
   }

   LinearVertexDistribution linear_vertex_distribution_;
};

TEST(LinearVertexDistributionTest0, AtLeastTwoVertices)
{
   LinearVertexDistribution lnd(0);

   EXPECT_EQ(2, lnd.number_of_vertices());
}

TEST_F(LinearVertexDistributionTest, ThrowForVertexOutOfRange)
{
   for (unsigned k = NumberOfVertices(); k < 2*NumberOfVertices(); ++k)
   {
      EXPECT_THROW(
         linear_vertex_distribution_.PositionNthVertex(k), OutOfRange);
   }
}

TEST_F(LinearVertexDistributionTest, NumberOfVertices)
{
   EXPECT_EQ(
      NumberOfVertices(), linear_vertex_distribution_.number_of_vertices());
}

TEST_F(LinearVertexDistributionTest, VertexPosition)
{
   for (unsigned k = NumberOfVertices()-1; k; --k)
   {
      EXPECT_FLOAT_EQ(
         linear_vertex_distribution_.PositionNthVertex(k),
         (.0+k)/(NumberOfVertices()-1.));
   }
}

// --- ExpVertexDistribution --------------------------------------
struct ExpVertexDistributionTest : testing::Test
{
   static constexpr std::size_t NumberOfVertices() {
      return 4;
   }
   static constexpr std::size_t kBase = 2;

   ExpVertexDistributionTest()
      : exp_vertex_distribution_(NumberOfVertices(), false, kBase)
   {
   }

   ExpVertexDistribution exp_vertex_distribution_;
};

TEST(ExpVertexDistributionTest0, InvalidBaseThrows)
{
   for (double base=.2; base <=1.; base+=.2)
   {
      EXPECT_THROW(
         ExpVertexDistribution(4, false, base), DomainError);
   }
}

TEST_F(ExpVertexDistributionTest, NumberOfVertices)
{
   EXPECT_EQ(
      NumberOfVertices(), exp_vertex_distribution_.number_of_vertices());
}

TEST_F(ExpVertexDistributionTest, ThrowForVertexOutOfRange)
{
   for (auto k=NumberOfVertices(); k < 2*NumberOfVertices(); ++k)
   {
      EXPECT_THROW(
         exp_vertex_distribution_.PositionNthVertex(k),
         OutOfRange);
   }
}

TEST_F(ExpVertexDistributionTest, VertexPosition)
{
   for (int i = NumberOfVertices() - 1; 0 < i; --i) {
      const double exponent = 1. + i - NumberOfVertices();
      const double power = std::pow(kBase, exponent);
      EXPECT_FLOAT_EQ(
         exp_vertex_distribution_.PositionNthVertex(i), power);
   }
   EXPECT_FLOAT_EQ(exp_vertex_distribution_.PositionNthVertex(0), .0);
}

TEST_F(ExpVertexDistributionTest, VertexPositionReverse)
{
   exp_vertex_distribution_.Reverse();

   for (std::size_t i = 0; NumberOfVertices() - 1 > i; ++i) {
      const double exponent = -static_cast<double>(i);
      const double power = 1. - std::pow(kBase, exponent);
      EXPECT_FLOAT_EQ(
         exp_vertex_distribution_.PositionNthVertex(i), power);
   }
   EXPECT_FLOAT_EQ(exp_vertex_distribution_.PositionNthVertex(
                      NumberOfVertices() - 1), 1.);
}

// --- GridPatch --------------------------------------------------

template <typename> struct GridPatchTestBase;

template <> struct GridPatchTestBase<CartesianSystem>
{
   using Point = CartesianSystem::Point;

   static constexpr std::size_t NumberOfVertexRows() {
      return 3;
   }
   static constexpr std::size_t NumberOfVertexColumns() {
      return 4;
   }

   Point LowerLeftCorner()
   {
      return Point(.5*metre, .5*metre);
   }

   Point UpperRightCorner()
   {
      return Point(2.*metre, 1.5*metre);
   }

   Length EdgeLenX(RowIndex, ColumnIndex)
   {
      using std::get;
      return (get<0>(UpperRightCorner())
              -get<0>(LowerLeftCorner()))/(-1.+NumberOfVertexColumns());
   }

   Length EdgeLenY(RowIndex, ColumnIndex)
   {
      using std::get;
      return (get<1>(UpperRightCorner())
              -get<1>(LowerLeftCorner()))/(-1.+NumberOfVertexRows());
   }

   Point EdgeMidPointX(RowIndex row_index, ColumnIndex column_index)
   {
      Length grid_extent_x
         = get<0>(UpperRightCorner()) - get<0>(LowerLeftCorner());
      Length grid_extent_y
         = get<1>(UpperRightCorner()) - get<1>(LowerLeftCorner());

      return Point(
         get<0>(LowerLeftCorner())
         + (.5+column_index.Get())/(NumberOfVertexColumns()-1.) * grid_extent_x,
         get<1>(LowerLeftCorner())
         + row_index.Get()/(NumberOfVertexRows()-1.) * grid_extent_y);
   }

   Point EdgeMidPointY(RowIndex row_index, ColumnIndex column_index)
   {
      Length grid_extent_x
         = get<0>(UpperRightCorner()) - get<0>(LowerLeftCorner());
      Length grid_extent_y
         = get<1>(UpperRightCorner()) - get<1>(LowerLeftCorner());

      return Point(
         get<0>(LowerLeftCorner())
         + (.0+column_index.Get())/(NumberOfVertexColumns()-1.) * grid_extent_x,
         get<1>(LowerLeftCorner())
         + (.5+row_index.Get())/(NumberOfVertexRows()-1.) * grid_extent_y);
   }

   Area CellSize(RowIndex row_index, ColumnIndex column_index)
   {
      return EdgeLenY(row_index, column_index)
         * EdgeLenX(row_index, column_index);
   }

   Length AverageEdgeLenX(
      RowIndex row_index, ColumnIndex column_index)
   {
      return EdgeLenX(row_index, column_index);
   }

   Length AverageEdgeLenY(
      RowIndex row_index, ColumnIndex column_index)
   {
      return EdgeLenY(row_index, column_index);
   }

   Length CrossingEdgeLengthX(
      RowIndex row_index, ColumnIndex column_index)
   {
      return EdgeLenY(row_index, column_index);
   }

   Length CrossingEdgeLengthY(
      RowIndex row_index, ColumnIndex column_index)
   {
      return EdgeLenX(row_index, column_index);
   }

   Point VertexPos(RowIndex row_index, ColumnIndex column_index)
   {
      using std::get;

      return Point(
         get<0>(LowerLeftCorner())+(get<0>(UpperRightCorner())
                                    -get<0>(LowerLeftCorner()))
         *(.0+column_index.Get())/(-1.+NumberOfVertexColumns()),
         get<1>(LowerLeftCorner())+(get<1>(UpperRightCorner())
                                    -get<1>(LowerLeftCorner()))
         *(.0+row_index.Get())/(-1.+NumberOfVertexRows()));
   }
};

template <> struct GridPatchTestBase<PolarSystem>
{
   using Point = PolarSystem::Point;

   static constexpr std::size_t NumberOfVertexRows() {
      return 3;
   }
   static constexpr std::size_t NumberOfVertexColumns() {
      return 4;
   }

   Point LowerLeftCorner()
   {
      return Point(.6*metre, .0*radians);
   }

   Point UpperRightCorner()
   {
      return Point(1.2*metre, kPi/2*radians);
   }

   Length Radius(ColumnIndex column_index)
   {
      return get<0>(LowerLeftCorner())
         + (get<0>(UpperRightCorner())-get<0>(LowerLeftCorner()))
         * (.0+column_index.Get())/(-1.+NumberOfVertexColumns());
   }

   Angle DeltaPhi()
   {
      return (
         get<1>(UpperRightCorner())
         -get<1>(LowerLeftCorner()))/(NumberOfVertexRows()-1.);
   }

   Length DeltaR()
   {
      return (
         get<0>(UpperRightCorner())
         -get<0>(LowerLeftCorner()))/(NumberOfVertexColumns()-1.);
   }

   Length EdgeLenX(RowIndex, ColumnIndex)
   {
      using std::get;
      return (get<0>(UpperRightCorner())
              -get<0>(LowerLeftCorner()))/(-1.+NumberOfVertexColumns());
   }

   Length EdgeLenY(RowIndex, ColumnIndex column_index)
   {
      using std::get;

      return Radius(column_index)*DeltaPhi()/radians;
   }

   Length AverageEdgeLenX(
      RowIndex row_index, ColumnIndex column_index)
   {
      return EdgeLenX(row_index, column_index);
   }

   Length AverageEdgeLenY(
      RowIndex row_index, ColumnIndex column_index)
   {
      return EdgeLenY(row_index, column_index);
   }

   Length CrossingEdgeLengthX(
      RowIndex row_index, ColumnIndex column_index)
   {
      const Length radius
         = (Radius(column_index)+Radius(column_index+1))/2.;

      return radius*DeltaPhi()/radians;
   }

   Length CrossingEdgeLengthY(
      RowIndex row_index, ColumnIndex column_index)
   {
      return EdgeLenX(row_index, column_index);
   }

   Area CellSize(RowIndex row_index, ColumnIndex column_index)
   {
      const Length r1 = (0 == column_index) ?
         Radius(column_index) - DeltaR()/2. :
         (Radius(column_index-1)+Radius(column_index))/2.;
      const Length r2 = (NumberOfVertexColumns()-1 == column_index) ?
         Radius(column_index) + DeltaR()/2. :
         (Radius(column_index)+Radius(column_index+1))/2.;

      return PolarSystem::SurfaceElement(
         Point(r1, .0*radians), Point(r2, DeltaPhi()));
   }

   Point EdgeMidPointX(RowIndex row_index, ColumnIndex column_index)
   {
      const Length r
         =(Radius(column_index)+Radius(column_index+1))/2.;

      const Angle phi = (.0+row_index.Get())*DeltaPhi();

      return Point(r, phi);
   }

   Point EdgeMidPointY(RowIndex row_index, ColumnIndex column_index)
   {
      const Length r = Radius(column_index);
      const Angle phi = (.5+row_index.Get())*DeltaPhi();

      return Point(r, phi);
   }

   Point VertexPos(RowIndex row_index, ColumnIndex column_index)
   {
      using std::get;

      return Point(
         get<0>(LowerLeftCorner())+(get<0>(UpperRightCorner())
                                    -get<0>(LowerLeftCorner()))
         *(.0+column_index.Get())/(-1.+NumberOfVertexColumns()),
         get<1>(LowerLeftCorner())+(get<1>(UpperRightCorner())
                                    -get<1>(LowerLeftCorner()))
         *(.0+row_index.Get())/(-1.+NumberOfVertexRows()));
   }
};

template <typename AnalysisAndCSys>
struct GridPatchTest
   : GridPatchTestBase<typename AnalysisAndCSys::CoordinateSystemType>,
   testing::Test
{
   using CoordinateSystem = typename AnalysisAndCSys::CoordinateSystemType;
   using Analysis = typename AnalysisAndCSys::AnalysisType;
   using LinearVertexDistributionPointer =
      std::shared_ptr<LinearVertexDistribution>;
   using GridPatchType = GridPatch<CoordinateSystem, Analysis>;
   using Time = typename Analysis::Time;

   GridPatchTest()
      : vertex_dist_x_(
         LinearVertexDistribution::New(this->NumberOfVertexColumns())),
        vertex_dist_y_(
           LinearVertexDistribution::New(this->NumberOfVertexRows())),
        grid_patch_(
           this->LowerLeftCorner(), this->UpperRightCorner(),
           vertex_dist_x_, vertex_dist_y_)
   {
   }

   LinearVertexDistributionPointer vertex_dist_x_;
   LinearVertexDistributionPointer vertex_dist_y_;
   GridPatch<CoordinateSystem, Analysis> grid_patch_;
};

TYPED_TEST_CASE(GridPatchTest, AnalysisAndCSysTypes);

TYPED_TEST(GridPatchTest, ThrowOnInvalidInitialization)
{
   using GridPatchType = typename GridPatchTest<TypeParam>::GridPatchType;

   EXPECT_THROW(GridPatchType(
                   this->UpperRightCorner(), this->LowerLeftCorner(),
                   this->vertex_dist_x_, this->vertex_dist_y_),
                InvalidArgument);
}

TYPED_TEST(GridPatchTest, MaxRowIndex)
{
   auto max_row_index = this->grid_patch_.MaxRowIndex();
   EXPECT_EQ(max_row_index.Get(), this->NumberOfVertexRows()-1);
}

TYPED_TEST(GridPatchTest, MaxColumnIndex)
{
   auto max_column_index = this->grid_patch_.MaxColumnIndex();
   EXPECT_EQ(max_column_index.Get(), this->NumberOfVertexColumns() -1);
}

TYPED_TEST(GridPatchTest, NumberOfVertexRows)
{
   EXPECT_EQ(this->NumberOfVertexRows(),
             this->grid_patch_.NumberOfVertexRows());
}

TYPED_TEST(GridPatchTest, NumberOfVertexColumns)
{
   EXPECT_EQ(this->NumberOfVertexColumns(),
             this->grid_patch_.NumberOfVertexColumns());
}

TYPED_TEST(GridPatchTest, EdgeLengthX)
{
   for (RowIndex row_index=RowIndex(0);
        row_index < this->NumberOfVertexRows(); ++row_index)
      for (ColumnIndex column_index=ColumnIndex(0);
           column_index < this->NumberOfVertexColumns()-1; ++column_index)
      {
         EXPECT_FLOAT_EQ(
            in_base_units(
               this->grid_patch_.EdgeLength(
                  row_index, column_index, Axis::kFirst)),
            in_base_units(this->EdgeLenX(row_index, column_index)))
            << "row index: " << row_index.Get()
            << " column index: " << column_index.Get();
      }
}

TYPED_TEST(GridPatchTest, EdgeLengthY)
{
   for (RowIndex row_index=RowIndex(0);
        row_index < this->NumberOfVertexRows()-1; ++row_index)
      for (ColumnIndex column_index=ColumnIndex(0);
           column_index < this->NumberOfVertexColumns(); ++column_index)
      {
         EXPECT_FLOAT_EQ(
            in_base_units(
               this->grid_patch_.EdgeLength(
                  row_index, column_index, Axis::kSecond)),
            in_base_units(this->EdgeLenY(row_index, column_index)))
            << "row index: " << row_index.Get()
            << " column index: " << column_index.Get();
      }
}

TYPED_TEST(GridPatchTest, AverageEdgeLength)
{
   for (RowIndex row_index=RowIndex(0);
        row_index < this->NumberOfVertexRows(); ++row_index)
      for (ColumnIndex column_index=ColumnIndex(0);
           column_index < this->NumberOfVertexColumns(); ++column_index)
      {
         EXPECT_FLOAT_EQ(
            in_base_units(
               this->grid_patch_.AverageEdgeLength(
                  row_index, column_index, Axis::kFirst)),
            in_base_units(this->AverageEdgeLenX(row_index, column_index)))
            << "row index: " << row_index.Get()
            << " column index: " << column_index.Get();

         EXPECT_FLOAT_EQ(
            in_base_units(
               this->grid_patch_.AverageEdgeLength(
                  row_index, column_index, Axis::kSecond)),
            in_base_units(this->AverageEdgeLenY(row_index, column_index)))
            << "row index: " << row_index.Get()
            << " column index: " << column_index.Get();
      }
}

TYPED_TEST(GridPatchTest, CrossingEdgeLength)
{
   for (RowIndex row_index=RowIndex(0);
        row_index < this->NumberOfVertexRows()-1; ++row_index)
      for (ColumnIndex column_index=ColumnIndex(0);
           column_index < this->NumberOfVertexColumns()-1; ++column_index)
      {
         EXPECT_FLOAT_EQ(
            in_base_units(
               this->grid_patch_.CrossingEdgeLength(
                  row_index, column_index, Axis::kFirst)),
            in_base_units(this->CrossingEdgeLengthX(row_index, column_index)));
         EXPECT_FLOAT_EQ(
            in_base_units(
               this->grid_patch_.CrossingEdgeLength(
                  row_index, column_index, Axis::kSecond)),
            in_base_units(this->CrossingEdgeLengthY(row_index, column_index)));
      }
}

TYPED_TEST(GridPatchTest, EdgeMidPointX)
{
   for (RowIndex row_index=RowIndex(0);
        row_index < this->NumberOfVertexRows()-1; ++row_index)
      for (ColumnIndex column_index=ColumnIndex(0);
           column_index < this->NumberOfVertexColumns()-1; ++column_index)
      {
         auto emx_grid = this->grid_patch_.EdgeMidPointLocal(
            row_index, column_index, Axis::kFirst);
         auto emx = this->EdgeMidPointX(row_index, column_index);
         EXPECT_FLOAT_EQ(
            in_base_units(get<0>(emx_grid)), in_base_units(get<0>(emx)));
         EXPECT_FLOAT_EQ(
            in_base_units(get<1>(emx_grid)), in_base_units(get<1>(emx)));
      }
}

TYPED_TEST(GridPatchTest, EdgeMidPointY)
{
   for (RowIndex row_index=RowIndex(0);
        row_index < this->NumberOfVertexRows()-1; ++row_index)
      for (ColumnIndex column_index=ColumnIndex(0);
           column_index < this->NumberOfVertexColumns()-1; ++column_index)
      {
         auto emy_grid = this->grid_patch_.EdgeMidPointLocal(
            row_index, column_index, Axis::kSecond);
         auto emy = this->EdgeMidPointY(row_index, column_index);
         EXPECT_FLOAT_EQ(
            in_base_units(get<0>(emy_grid)), in_base_units(get<0>(emy)));
         EXPECT_FLOAT_EQ(
            in_base_units(get<1>(emy_grid)), in_base_units(get<1>(emy)));
      }
}

TYPED_TEST(GridPatchTest, CellSize)
{
   for (RowIndex row_index=RowIndex(0);
        row_index < this->NumberOfVertexRows(); ++row_index)
      for (ColumnIndex column_index=ColumnIndex(0);
           column_index < this->NumberOfVertexColumns(); ++column_index)
      {
         EXPECT_FLOAT_EQ(
            in_base_units(
               this->grid_patch_.CellSize(row_index, column_index)),
            in_base_units(
               this->CellSize(row_index, column_index)));
      }
}

TYPED_TEST(GridPatchTest, VertexPositionLocal)
{
   using std::get;

   for (RowIndex row_index=RowIndex(0);
        row_index < this->NumberOfVertexRows(); ++row_index)
      for (ColumnIndex column_index=ColumnIndex(0);
           column_index < this->NumberOfVertexColumns(); ++column_index)
      {
         auto vp1 = this->grid_patch_.VertexPositionLocal(
            row_index, column_index);
         auto vp2 = this->VertexPos(row_index, column_index);

         EXPECT_FLOAT_EQ(
            in_base_units(get<0>(vp1)), in_base_units(get<0>(vp2)));
         EXPECT_FLOAT_EQ(
            in_base_units(get<1>(vp1)), in_base_units(get<1>(vp2)));
      }
}

TYPED_TEST(GridPatchTest, GetSetDeleteMaterial)
{
   using Point = typename TypeParam::CoordinateSystemType::Point;
   using Rect = Rectangle<typename TypeParam::CoordinateSystemType>;

   const Point grid_midpoint = ManhattenMidpoint(
      this->LowerLeftCorner(), this->UpperRightCorner());

   this->grid_patch_.DeleteMaterial(grid_midpoint);

   EXPECT_FLOAT_EQ(
      in_base_units(
         get<0>(
            this->grid_patch_.GetMaterial(grid_midpoint)->EvaluateReluctivity(
               1.*tesla, Axis::kFirst))),
      in_base_units(1./mu_0));

   this->grid_patch_.SetMaterial(
      Rect::New(this->LowerLeftCorner(), this->UpperRightCorner()),
      Material::New(
         _isotropic_reluctivity = ConstantReluctivity::New(2.)));

   EXPECT_FLOAT_EQ(
      in_base_units(
         get<0>(
            this->grid_patch_.GetMaterial(grid_midpoint)->EvaluateReluctivity(
               1.*tesla, Axis::kFirst))),
      in_base_units(.5/mu_0));

   this->grid_patch_.DeleteMaterial(grid_midpoint);

   EXPECT_FLOAT_EQ(
      in_base_units(
         get<0>(
            this->grid_patch_.GetMaterial(grid_midpoint)->EvaluateReluctivity(
               1.*tesla, Axis::kFirst))),
      in_base_units(1./mu_0));
}

TYPED_TEST(GridPatchTest, GetSetDeleteLoad)
{
   using Point = typename TypeParam::CoordinateSystemType::Point;
   using Rect = Rectangle<typename TypeParam::CoordinateSystemType>;
   using Time = typename TypeParam::AnalysisType::Time;

   const Point grid_midpoint = ManhattenMidpoint(
      this->LowerLeftCorner(), this->UpperRightCorner());

   this->grid_patch_.DeleteLoad(grid_midpoint);

   EXPECT_FLOAT_EQ(
      abs(in_base_units(
             (* (this->grid_patch_.GetLoad(grid_midpoint)))(
                Time::from_value(.0)))), .0);

   this->grid_patch_.SetLoad(
      Rect::New(this->LowerLeftCorner(), this->UpperRightCorner()),
      1.*ampere_per_metre_squared);

   EXPECT_FLOAT_EQ(
      abs(in_base_units(
             (* (this->grid_patch_.GetLoad(grid_midpoint)))(
                Time::from_value(.0)))), 1.);

   this->grid_patch_.DeleteLoad(grid_midpoint);

   EXPECT_FLOAT_EQ(
      abs(in_base_units(
             (* (this->grid_patch_.GetLoad(grid_midpoint)))(
                Time::from_value(.0)))), .0);
}

TYPED_TEST(GridPatchTest, RotateShift)
{
   this->grid_patch_.Rotate(kPi/2.*radians);
   const CartesianSystem::Vector shift_vector(-1. * metre, 1. * metre);
   this->grid_patch_.Shift(shift_vector);

   auto p = this->grid_patch_.VertexPositionGlobal(RowIndex(0), ColumnIndex(0));
   auto llc_global = TypeParam::CoordinateSystemType::ToCartesic(this->LowerLeftCorner());
   EXPECT_FLOAT_EQ(
      in_base_units(get<0>(p)),
      in_base_units(-get<1>(llc_global)+get<0>(shift_vector)));
   EXPECT_FLOAT_EQ(
      in_base_units(get<1>(p)),
      in_base_units(get<0>(llc_global)+get<1>(shift_vector)));
}

// --- Simple coupling -----------------------------------------------

template <typename Analysis>
struct SimpleCouplingTest : testing::Test
{
   using GridPatchCartesian = GridPatch<CartesianSystem, Analysis>;
   using GPCPointer = std::shared_ptr<GridPatchCartesian>;
   using LVDPointer = std::shared_ptr<LinearVertexDistribution>;
   using GAPointer = std::shared_ptr<GridAssembly>;
   using Point = CartesianSystem::Point;

   static Point P1()
   {
      return Point(.0 * metre, .0 * metre);
   }
   static Point P2()
   {
      return Point(1. * metre, 1. * metre);
   }
   static Point P3()
   {
      return Point(1. * metre, .0 * metre);
   }
   static Point P4()
   {
      return Point(1.5 * metre, 1. * metre);
   }

   SimpleCouplingTest()
      : vertex_dist_(LinearVertexDistribution::New(4)),
        grid_patches_({
              GridPatchCartesian::New(
                 P1(), P2(), vertex_dist_, vertex_dist_),
                 GridPatchCartesian::New(
                    P3(), P4(), vertex_dist_, vertex_dist_)}),
        grid_assembly_(GridAssembly::New(
                          &(grid_patches_[0]), &(grid_patches_[2])))

   {
      grid_patches_[0]->east_boundary().SetBoundaryCoupling(
         grid_patches_[1]->west_boundary());
      grid_assembly_->SetFirstVertexIndex();
   }

   LVDPointer vertex_dist_;
   GPCPointer grid_patches_[2];
   GAPointer grid_assembly_;
};

TYPED_TEST_CASE(SimpleCouplingTest, AnalysisTypes);

TYPED_TEST(SimpleCouplingTest, AverageEdgeLength)
{
   EXPECT_FLOAT_EQ(
      in_base_units(
         this->grid_patches_[0]->AverageEdgeLength(
            RowIndex(1), ColumnIndex(3), Axis::kFirst)), .25);
}

TYPED_TEST(SimpleCouplingTest, CellSize)
{
   EXPECT_FLOAT_EQ(
      in_base_units(
         this->grid_patches_[0]->CellSize(
            RowIndex(1), ColumnIndex(3))), 1./12);
}

TYPED_TEST(SimpleCouplingTest, CrossingEdgeLength)
{
   EXPECT_FLOAT_EQ(
      in_base_units(
         this->grid_patches_[0]->CrossingEdgeLength(
            RowIndex(1), ColumnIndex(3), Axis::kSecond)), .25);
}

TYPED_TEST(SimpleCouplingTest, GlobalIndexAtBoundary)
{
   EXPECT_EQ(
      this->grid_patches_[0]->GlobalIndex(RowIndex(1), ColumnIndex(3)), 20);
}

// --- Circular coupling ---------------------------------------------
template <typename Analysis>
struct CircularCouplingTest : testing::Test
{
   using GridPatchCartesian = GridPatch<CartesianSystem, Analysis>;
   using GPCPointer = std::shared_ptr<GridPatchCartesian>;
   using LVDPointer = std::shared_ptr<LinearVertexDistribution>;
   using Point = CartesianSystem::Point;

   static Point P1() {return Point(.0 * metre, .0 * metre);}
   static Point P2() {return Point(1. * metre, 1. * metre);}
   static Point P3() {return Point(1. * metre, .0 * metre);}
   static Point P4() {return Point(2. * metre, 1. * metre);}
   static Point P5() {return Point(1. * metre, 1. * metre);}
   static Point P6() {return Point(2. * metre, 2. * metre);}
   static Point P7() {return Point(0. * metre, 1. * metre);}
   static Point P8() {return Point(1. * metre, 2. * metre);}

   CircularCouplingTest()
      : vertex_dist_x_(LinearVertexDistribution::New(4)),
        vertex_dist_y_(LinearVertexDistribution::New(3)),
        grid_patches_({
              GridPatchCartesian::New(
                 P1(), P2(), vertex_dist_x_, vertex_dist_y_),
                 GridPatchCartesian::New(
                    P3(), P4(), vertex_dist_x_, vertex_dist_y_),
                 GridPatchCartesian::New(
                    P5(), P6(), vertex_dist_x_, vertex_dist_y_),
                 GridPatchCartesian::New(
                    P7(), P8(), vertex_dist_x_, vertex_dist_y_)
                 }),
        grid_assembly_(grid_patches_, grid_patches_+4)
   {
      grid_patches_[0]->east_boundary().SetBoundaryCoupling(
         grid_patches_[1]->west_boundary());
      grid_patches_[1]->north_boundary().SetBoundaryCoupling(
         grid_patches_[2]->south_boundary());
      grid_patches_[2]->west_boundary().SetBoundaryCoupling(
         grid_patches_[3]->east_boundary());
      grid_patches_[3]->south_boundary().SetBoundaryCoupling(
         grid_patches_[0]->north_boundary());

      grid_assembly_.SetFirstVertexIndex();
   }

   LVDPointer vertex_dist_x_;
   LVDPointer vertex_dist_y_;
   GPCPointer grid_patches_[4];
   GridAssembly grid_assembly_;
};

TYPED_TEST_CASE(CircularCouplingTest, AnalysisTypes);

TYPED_TEST(CircularCouplingTest, GlobalIndexAtCenterPoint)
{
   EXPECT_EQ(
      this->grid_patches_[0]->GlobalIndex(RowIndex(2), ColumnIndex(3)), 11);
}

// --- GraphProperties -------------------------------------------
template <typename Analysis>
struct GraphPropertiesTest : testing::Test
{
   using Time = typename Analysis::Time;

   std::size_t kNumberOfEquations = 1234;

   GraphPropertiesTest()
      : graph_properties_(kNumberOfEquations, Time::from_value(1.))
   {
   }

   GraphProperties<Analysis> graph_properties_;
};

TYPED_TEST_CASE(GraphPropertiesTest, AnalysisTypes);

TYPED_TEST(GraphPropertiesTest, NumberOfEquations) {
   EXPECT_EQ(this->graph_properties_.number_of_equations(),
             this->kNumberOfEquations);
   this->graph_properties_.set_number_of_equations(
      2*(this->kNumberOfEquations));
   EXPECT_EQ(this->graph_properties_.number_of_equations(),
             2*(this->kNumberOfEquations));
}

TYPED_TEST(GraphPropertiesTest, TimeStep) {
   using Time = typename TypeParam::Time;
   EXPECT_FLOAT_EQ(in_base_units(this->graph_properties_.time_step()), 1.);
   this->graph_properties_.set_time_step(Time::from_value(2.));
   EXPECT_FLOAT_EQ(in_base_units(this->graph_properties_.time_step()), 2.);
}

TYPED_TEST(GraphPropertiesTest, TimeStepping) {

   for (int i = 0; i < 10; ++i)
   {
      EXPECT_FLOAT_EQ(in_base_units(this->graph_properties_.current_time()), i);
      this->graph_properties_.NextTimeStep();
   }
}

// --- fi::StraightGridEdge ------------------------------------------
template <typename Analysis>
struct StraightGridEdgeTest : testing::Test
{
   using MagneticFluxDensity = typename Analysis::MagneticFluxDensity;

   StraightGridEdgeTest()
      : constant_reluctivity_(ConstantReluctivity::New(1.)),
        straight_grid_edge_(.1*metre, .5, constant_reluctivity_)
   {
   }

   ConstReluctivityPointer constant_reluctivity_;
   fi::StraightGridEdge<Analysis> straight_grid_edge_;
};

TYPED_TEST_CASE(StraightGridEdgeTest, AnalysisTypes);

TYPED_TEST(StraightGridEdgeTest, GetSetFluxDensity) {
   this->straight_grid_edge_.set_flux_density(
      TypeParam::MagneticFluxDensity::from_value(1.));

   EXPECT_FLOAT_EQ(1., in_base_units(
                      abs(this->straight_grid_edge_.flux_density())));
}

TYPED_TEST(StraightGridEdgeTest, MagneticVoltage) {
   this->straight_grid_edge_.set_flux_density(
      TypeParam::MagneticFluxDensity::from_value(1.));

   EXPECT_FLOAT_EQ(
      in_base_units(1.*tesla/mu_0*.1*metre*.5),
      in_base_units(abs(this->straight_grid_edge_.magnetic_voltage())));
}

// --- fi::CircularGridEdge ------------------------------------------
template <typename Analysis>
struct CircularGridEdgeTest : testing::Test
{
   CircularGridEdgeTest()
      : constant_reluctivity_(ConstantReluctivity::New(1.)),
        circular_grid_edge_(.1*metre, .5, .4*metre, constant_reluctivity_)
   {
   }

   ConstReluctivityPointer constant_reluctivity_;
   fi::CircularGridEdge<Analysis> circular_grid_edge_;
};

TYPED_TEST_CASE(CircularGridEdgeTest, AnalysisTypes);

TYPED_TEST(CircularGridEdgeTest, GetSetFluxDensity) {
   using MagneticFluxDensity = typename TypeParam::MagneticFluxDensity;

   this->circular_grid_edge_.set_flux_density(
      MagneticFluxDensity::from_value(1.));

   EXPECT_FLOAT_EQ(1., in_base_units(
                      abs(this->circular_grid_edge_.flux_density())));
}

TYPED_TEST(CircularGridEdgeTest, MagneticVoltage) {
   using MagneticFluxDensity = typename TypeParam::MagneticFluxDensity;
   this->circular_grid_edge_.set_flux_density(
      MagneticFluxDensity::from_value(1.));

   EXPECT_FLOAT_EQ(
      in_base_units(1.*tesla/mu_0*.1*metre*.5),
      in_base_units(abs(this->circular_grid_edge_.magnetic_voltage())));
}

// --- StandardVertex ------------------------------------------------
template <typename Analysis>
struct StandardVertexAttributes
{
   using Point = CartesianSystem::Point;
   using CurrentDensity = typename Analysis::CurrentDensity;

   static Point Position() {
      return Point(.5*metre, -.4*metre);
   }

   static Area Surface() {
      return .42*metre*metre;
   }

   static CurrentDensity CD() {
      return 5.6 * ampere_per_metre_squared;
   }
};

template <typename Analysis>
struct StandardVertexTest
   : StandardVertexAttributes<Analysis>, testing::Test
{
   using Point = CartesianSystem::Point;
   using RhsVector = typename Analysis::RhsVector;

   StandardVertexTest()
      : load_definition_(LoadDefinition<Analysis>::New(this->CD())),
        standard_vertex_(this->Position(), load_definition_, this->Surface())
   {
   }

   ConstLoadDefinitionPointer<Analysis> load_definition_;
   fi::StandardVertex<Analysis> standard_vertex_;
};

TYPED_TEST_CASE(StandardVertexTest, AnalysisTypes);

TYPED_TEST(StandardVertexTest, SetGetEquationNumber) {
   this->standard_vertex_.set_equation_number(12345);
   EXPECT_EQ(this->standard_vertex_.equation_number(), 12345);
}

TYPED_TEST(StandardVertexTest, GetPosition) {
   auto position = this->standard_vertex_.position();

   EXPECT_FLOAT_EQ(in_base_units(get<0>(position)),
                   in_base_units(get<0>(this->Position())));

   EXPECT_FLOAT_EQ(in_base_units(get<1>(position)),
                   in_base_units(get<1>(this->Position())));
}

TYPED_TEST(StandardVertexTest, SetPosition) {
   using Point = CartesianSystem::Point;

   this->standard_vertex_.set_position(
      Point(get<1>(this->Position()), get<0>(this->Position())));

   auto position = this->standard_vertex_.position();

   EXPECT_FLOAT_EQ(in_base_units(get<0>(position)),
                   in_base_units(get<1>(this->Position())));

   EXPECT_FLOAT_EQ(in_base_units(get<1>(position)),
                   in_base_units(get<0>(this->Position())));
}

TYPED_TEST(StandardVertexTest, SetGetPotential) {
   typename TypeParam::RhsVector rhs_vector(3);
   rhs_vector.fill(1.1);

   this->standard_vertex_.set_equation_number(2);
   this->standard_vertex_.set_potential(rhs_vector);

   EXPECT_FLOAT_EQ(in_base_units(
                      abs(this->standard_vertex_.potential())), 1.1);
}

// --- NeumannVertex -------------------------------------------------
template <typename Analysis>
struct NeumannVertexTest
   : StandardVertexAttributes<Analysis>, testing::Test
{
   using LoadDefPtr = ConstLoadDefinitionPointer<Analysis>;
   using Point = CartesianSystem::Point;
   using MagneticFluxDensity = typename Analysis::MagneticFluxDensity;
   using RhsVector = typename Analysis::RhsVector;

   static MagneticReluctivity ReluctOutside() {
      return 10./mu_0;
   }

   static MagneticFluxDensity FluxDensOut() {
      return .8*tesla;
   }

   static Length GridSpacing() {
      return .4e-2*metre;
   }

   NeumannVertexTest()
      : load_definition_(LoadDefinition<Analysis>::New(this->CD())),
        neumann_vertex_(
           this->Position(), this->ReluctOutside(),
           this->FluxDensOut(), this->GridSpacing(),
           Parallelity::kParallel, load_definition_, this->Surface())
   {
   }

   LoadDefPtr load_definition_;
   fi::NeumannVertex<Analysis> neumann_vertex_;
};

TYPED_TEST_CASE(NeumannVertexTest, AnalysisTypes);

TYPED_TEST(NeumannVertexTest, SetGetEquationNumber) {
   this->neumann_vertex_.set_equation_number(12345);
   EXPECT_EQ(this->neumann_vertex_.equation_number(), 12345);
}

TYPED_TEST(NeumannVertexTest, GetPosition) {
   auto position = this->neumann_vertex_.position();

   EXPECT_FLOAT_EQ(in_base_units(get<0>(position)),
                   in_base_units(get<0>(this->Position())));

   EXPECT_FLOAT_EQ(in_base_units(get<1>(position)),
                   in_base_units(get<1>(this->Position())));
}

TYPED_TEST(NeumannVertexTest, SetPosition) {
   using Point = CartesianSystem::Point;

   this->neumann_vertex_.set_position(
      Point(get<1>(this->Position()), get<0>(this->Position())));

   auto position = this->neumann_vertex_.position();

   EXPECT_FLOAT_EQ(in_base_units(get<0>(position)),
                   in_base_units(get<1>(this->Position())));

   EXPECT_FLOAT_EQ(in_base_units(get<1>(position)),
                   in_base_units(get<0>(this->Position())));
}

TYPED_TEST(NeumannVertexTest, SetGetPotential) {
   typename TypeParam::RhsVector rhs_vector(3);
   rhs_vector.fill(1.1);

   this->neumann_vertex_.set_equation_number(2);
   this->neumann_vertex_.set_potential(rhs_vector);

   EXPECT_FLOAT_EQ(in_base_units(abs(this->neumann_vertex_.potential())), 1.1);
}

// --- DirichletVertex -----------------------------------------------
template <typename Analysis>
struct DirichletVertexTest : testing::Test
{
   using MagneticVectorPotential = typename Analysis::MagneticVectorPotential;
   using Point = CartesianSystem::Point;

   static Point Position() {
      return Point(.5*metre, -.4*metre);
   }

   static MagneticVectorPotential Potential() {
      return MagneticVectorPotential::from_value(1e-3);
   }

   DirichletVertexTest()
      : dirichlet_vertex_(Position(), Potential())
   {
   }

   fi::DirichletVertex<Analysis> dirichlet_vertex_;
};

TYPED_TEST_CASE(DirichletVertexTest, AnalysisTypes);

TYPED_TEST(DirichletVertexTest, GetPosition) {
   auto position = this->dirichlet_vertex_.position();

   EXPECT_FLOAT_EQ(in_base_units(get<0>(position)),
                   in_base_units(get<0>(this->Position())));

   EXPECT_FLOAT_EQ(in_base_units(get<1>(position)),
                   in_base_units(get<1>(this->Position())));
}

TYPED_TEST(DirichletVertexTest, SetPosition) {
   using Point = CartesianSystem::Point;

   this->dirichlet_vertex_.set_position(
      Point(get<1>(this->Position()), get<0>(this->Position())));

   auto position = this->dirichlet_vertex_.position();

   EXPECT_FLOAT_EQ(in_base_units(get<0>(position)),
                   in_base_units(get<1>(this->Position())));

   EXPECT_FLOAT_EQ(in_base_units(get<1>(position)),
                   in_base_units(get<0>(this->Position())));
}

TYPED_TEST(DirichletVertexTest, Potential) {
   EXPECT_FLOAT_EQ(in_base_units(abs(this->dirichlet_vertex_.potential())),
                   in_base_units(abs(this->Potential())));
}

// --- CompileSimpleGrid ---------------------------------------------
template <typename> struct CompileSimpleGridBase;

template <> struct CompileSimpleGridBase<CartesianSystem> {
   using Point = CartesianSystem::Point;

   static Point LowerLeftCorner() {
      return Point(.0*metre, .0*metre);
   }

   static Point UpperRightCorner() {
      return Point(1.1*metre, 1.1*metre);
   }
};

template <> struct CompileSimpleGridBase<PolarSystem> {
   using Point = PolarSystem::Point;

   static Point LowerLeftCorner() {
      return Point(.05*metre, .0*radians);
   }

   static Point UpperRightCorner() {
      return Point(1.1*metre, 1.*radians);
   }
};

template <typename AnalysisAndCS>
struct CompileSimpleGridTest
   : CompileSimpleGridBase<typename AnalysisAndCS::CoordinateSystemType>,
   testing::Test
{
   using Analysis = typename AnalysisAndCS::AnalysisType;
   using CoordinateSystem = typename AnalysisAndCS::CoordinateSystemType;
   using GridType = GridPatch<CoordinateSystem, Analysis>;
   using GridPatchPointer = std::shared_ptr<GridType>;
   using LVDPointer = std::shared_ptr<LinearVertexDistribution>;
   using Time = typename Analysis::Time;

   CompileSimpleGridTest()
      : vertex_dist_(LinearVertexDistribution::New(4)),
        grid_patch_(
           GridType::New(this->LowerLeftCorner(), this->UpperRightCorner(),
                         vertex_dist_, vertex_dist_)),
        grid_compiler_(grid_patch_)
   {
   }

   LVDPointer vertex_dist_;
   GridPatchPointer grid_patch_;
   fi::GridCompiler<Analysis> grid_compiler_;
};

TYPED_TEST_CASE(CompileSimpleGridTest, AnalysisAndCSysTypes);

TYPED_TEST(CompileSimpleGridTest, NumberOfVertices) {
   using Time = typename TypeParam::AnalysisType::Time;

   auto grid_graph = this->grid_compiler_(Time::from_value(1.1));
   // WriteToGraphviz(grid_graph, "Simple_Grid.dot");
   EXPECT_EQ(16, num_vertices(*grid_graph));
}

// --- GridCompiler --------------------------------------------------
template <typename Analysis>
struct GridCompilerTest : SimpleCouplingTest<Analysis>
{
   using Time = typename Analysis::Time;

   static Time TimeStep() {
      return Time::from_value(1e-3);
   }

   GridCompilerTest() : grid_compiler_(this->grid_assembly_)
   {
   }

   fi::GridCompiler<Analysis> grid_compiler_;
};

TYPED_TEST_CASE(GridCompilerTest, AnalysisTypes);

TYPED_TEST(GridCompilerTest, NumberOfVertices)
{
   auto grid_graph = this->grid_compiler_(this->TimeStep());
   // WriteToGraphviz(grid_graph, "Simple_Coupling.dot");
   EXPECT_EQ(num_vertices(*grid_graph), 28);
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
