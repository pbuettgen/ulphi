// -*- coding: utf-8 -*-
/*
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and library for
 * computing electromagnetic fields.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "config.h"
#include "coordinate_systems.hpp"
#include "types.hpp"
#include "quantity_type_conversion.hpp"
#include "scalar_reluctivity.hpp"
#include "tuple_type_conversion.hpp"

//! Initialize type conversion for quantities
/*!
 *  Initialize type conversion for quantities with float type
 *  selection.
 *
 *  \tparam Float is the desired float type.
 */
template <typename Float>
void QuantityTypeConversion() {
   using namespace boost::units;
   using namespace PACKAGE_NS::python;

   QuantityRegisterConversion<si::current, Float>();
   QuantityRegisterConversion<si::current_density, Float>();
   QuantityRegisterConversion<si::magnetic_field_intensity, Float>();
   QuantityRegisterConversion<si::magnetic_flux_density, Float>();
   QuantityRegisterConversion<si::magnetic_vector_potential, Float>();
}

template <typename CoordinateSystem>
inline void RegisterPositionVectorConversion()
{
   PACKAGE_NS::python::TupleRegisterConversion<
      typename CoordinateSystem::FirstCoordinate,
      typename CoordinateSystem::SecondCoordinate>();
}

BOOST_PYTHON_MODULE(MODULE_NAME) {
   using namespace PACKAGE_NS;
   using namespace PACKAGE_NS::python;
   using namespace boost::units;

   QuantityTypeConversion<double>();
   QuantityTypeConversion<std::complex<double> >();

   QuantityRegisterConversion((ElectricalConductivity*) nullptr);
   QuantityRegisterConversion((Frequency*) nullptr);
   QuantityRegisterConversion((Length*) nullptr);
   QuantityRegisterConversion((MagneticReluctivity*) nullptr);
   QuantityRegisterConversion((MagneticReluctivityDerived*) nullptr);
   QuantityRegisterConversion((Permeability*) nullptr);
   QuantityRegisterConversion((Angle*) nullptr);
   QuantityRegisterConversion((Time*) nullptr);

   TupleRegisterConversion<
      ConstReluctivityPointer, ConstReluctivityPointer>();
   TupleRegisterConversion<
      MagneticReluctivity, MagneticReluctivityDerived>();
   RegisterPositionVectorConversion<CartesianSystem>();
   RegisterPositionVectorConversion<PolarSystem>();
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
