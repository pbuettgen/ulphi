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
 *  \brief type definitions
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef TYPES_HPP_nkr54TI4_
#define TYPES_HPP_nkr54TI4_

#include <tuple>
#include <vector>

#include "config.h"

#include <boost/units/systems/si.hpp>

#include <Eigen/Sparse>

#include "boost_units_extra.hpp"

namespace PACKAGE_NS {
   //! Angle
   using Angle = boost::units::quantity<boost::units::si::plane_angle, double>;

   //! Conductivity type
   using ElectricalConductivity =
      boost::units::quantity<boost::units::si::conductivity, double>;

   //! Frequency type
   using Frequency =
      boost::units::quantity<boost::units::si::frequency, double>;

   //! Permeability type
   using Permeability =
      boost::units::quantity<boost::units::si::permeability, double>;

   //! Magnetic reluctivity type
   using MagneticReluctivity =
      boost::units::quantity<boost::units::si::magnetic_reluctivity, double>;

   /*!
    *  Type for the magnetic reluctivity derived with respect to the
    *  magnetic flux density.
    */
   using MagneticReluctivityDerived = boost::units::quantity<
      boost::units::si::magnetic_reluctivity_derived, double>;

   //! Reluctivity and its derivation
   using ReluctivityReluctivityDerivedPair =
      std::tuple<MagneticReluctivity, MagneticReluctivityDerived>;

   //! Magnetic field intensity type (real value type)
   using MagneticFieldIntensityReal = boost::units::quantity<
      boost::units::si::magnetic_field_intensity, double>;

   //! Magnetic flux density type (real value type)
   using MagneticFluxDensityReal = boost::units::quantity<
      boost::units::si::magnetic_flux_density, double>;

   //! Magnetic polarization type
   using MagneticPolarizationReal = MagneticFluxDensityReal;

   //! Permittivity
   using Permittivity = boost::units::quantity<
      boost::units::si::permittivity, double>;

   //! Thermal conductivity
   using ThermalConductivity = boost::units::quantity<
      boost::units::si::thermal_conductivity, double>;

   //! Time type
   using Time =
      boost::units::quantity<boost::units::si::time, double>;

   //! Vertex index type
   using VertexIndex = std::size_t;

   //! Equation number type
   using EquationNumber = std::size_t;

   //! Length type
   using Length = boost::units::quantity<boost::units::si::length, double>;

   //! Area type
   using Area = boost::units::quantity<boost::units::si::area, double>;
}

#endif // TYPES_HPP_nkr54TI4_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
