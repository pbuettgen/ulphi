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
#include <tuple>
#include <random>
#include <iostream>

#include "config.h"

#include <boost/units/io.hpp>

#include <gtest/gtest.h>

#include "analysis.hpp"
#include "types.hpp"
#include "python/quantity_type_conversion.hpp"
#include "python/tuple_type_conversion.hpp"

using namespace PACKAGE_NS;
using namespace PACKAGE_NS::python;
namespace bp = boost::python;

#define ANALYSIS_DEPENDEND_TYPES(ANALYSIS_TYPE)             \
   typename ANALYSIS_TYPE<>::Current,                       \
      typename ANALYSIS_TYPE<>::CurrentDensity,             \
      typename ANALYSIS_TYPE<>::MagneticFieldIntensity,     \
      typename ANALYSIS_TYPE<>::MagneticFluxDensity,        \
      typename ANALYSIS_TYPE<>::MagneticVectorPotential,    \
      typename ANALYSIS_TYPE<>::Time

using QuantityTypes = testing::Types<
   ElectricalConductivity,
   Length,
   Permeability,
   MagneticReluctivity,
   MagneticReluctivityDerived,
   ANALYSIS_DEPENDEND_TYPES(MagnetoHarmonicAnalysis),
   ANALYSIS_DEPENDEND_TYPES(TimeTransientAnalysis)
   >;

#undef ANALYSIS_DEPENDEND_TYPES

struct PythonInitializeFinalize {
   PythonInitializeFinalize() {
      Py_Initialize();
   }
   ~PythonInitializeFinalize() {
      Py_Finalize();
   }
} instance;

double MakeRandom(double*) {
   static std::default_random_engine rand_source;
   static std::uniform_real_distribution<> urd(
      std::numeric_limits<float>::min(),
      std::numeric_limits<float>::max());

   return urd(rand_source);
}

std::complex<double> MakeRandom(std::complex<double>*) {
   return std::complex<double>(
      MakeRandom((double*)nullptr), MakeRandom((double*)nullptr));
}

// --- Python unit type conversion ----------------------------------

template <typename Quantity>
struct PythonUnitTypeConversionTest : testing::Test
{
   PythonUnitTypeConversionTest()
      : some_quantity_(
         Quantity::from_value(
            MakeRandom(
               static_cast<typename Quantity::value_type*>(nullptr))))
   {}

   Quantity some_quantity_;
};

TYPED_TEST_CASE(PythonUnitTypeConversionTest, QuantityTypes);

TYPED_TEST(PythonUnitTypeConversionTest, GeneralTest)
{
   using Unit = typename TypeParam::unit_type;
   using Float = typename TypeParam::value_type;

   bp::object python_quantity(
      bp::handle<>(
         QuantityToPythonConversion<Unit, Float>::convert(
            this->some_quantity_)));

   const auto is_convertible
      = QuantityIsConvertible<Unit, Float>(python_quantity.ptr());

   EXPECT_EQ(is_convertible, python_quantity.ptr());

   auto cpp_quantity = python::detail::ConvertFromPythonPint<Unit, Float>(
      python_quantity.ptr());

   EXPECT_EQ(cpp_quantity, this->some_quantity_);
}

// --- Tuple type conversion -----------------------------------------
struct TupleTypeConversionTest : testing::Test
{
   using Tuple = std::tuple<char, unsigned>;

   TupleTypeConversionTest()
      : t_(std::make_tuple('A', 1111))
   {}

   Tuple t_;
};

TEST_F(TupleTypeConversionTest, GeneralTest)
{
   bp::object python_tuple(
      bp::handle<>(
         TupleToPythonConversion<char, unsigned>::convert(
            this->t_)));

   const auto is_convertible
      =TupleIsConvertible<char, unsigned>(python_tuple.ptr());

   EXPECT_EQ(is_convertible, python_tuple.ptr());

   python::detail::TuplePythonToCppWrapper<char, unsigned> wrapper(
      python_tuple.ptr());

   EXPECT_EQ(this->t_, wrapper());
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
