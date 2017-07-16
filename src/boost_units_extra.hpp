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

/*! \file
 *
 *  \brief extra units for the boost units library
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 *
 *  Definition of extra units for the boost units library.  See the <a
 *  href="http://localhost/cgi-bin/dwww/usr/share/doc/libboost1.49-doc/HTML/doc/html/boost_units.html">boost
 *  units library's documentation</a> for more details.
 */

#ifndef BOOST_UNITS_EXTRA_HPP_FH1NGMhc_
#define BOOST_UNITS_EXTRA_HPP_FH1NGMhc_

#include <cmath>
#include <complex>

#include <boost/functional/hash.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/derived_dimension.hpp>
#include <boost/units/dimensionless_type.hpp>
#include <boost/units/physical_dimensions/current.hpp>
#include <boost/units/physical_dimensions/length.hpp>
#include <boost/units/physical_dimensions/mass.hpp>
#include <boost/units/physical_dimensions/time.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/unit.hpp>

namespace boost {
   namespace units {
      //! Compute a complex quantity's magnitude
      /*!
       *  \tparam Unit is the unit type
       *
       *  \tparam Y is some numeric type
       *
       *  \param source is the source quantity
       *
       *  \returns the argument's magnitude
       */
      template<typename Unit, typename Y>
      quantity<Unit, Y> abs(const quantity<Unit, std::complex<Y> >& source)
      {
         return quantity<Unit, Y>::from_value(std::abs(source.value()));
      }

      //! Compute a quantity's hash value
      /*!
       *  The computed hash value is computed from the quantity's
       *  value.
       *
       *  \param source is the quantity for which the hash value shall
       *  be computed.
       *
       *  \returns a hash value for the quantity passed in as paramter.
       */
      template<typename Unit, typename Y>
      std::size_t hash_value(const quantity<Unit, Y>& source)
      {
         boost::hash<Y> hasher;
         return hasher(source.value());
      }

      template <typename Unit, typename X, typename Y>
      bool FloatEqual(const quantity<Unit, X>& q1,
                      const quantity<Unit, Y>& q2)
      {
         return FloatEqual(in_base_units(q1), in_base_units(q2));
      }

      //! current density dimension
      using current_density_dimension =
         derived_dimension<length_base_dimension, -2,
                           current_base_dimension, 1>::type;

      //! magnetic vector potential dimension
      using magnetic_vector_potential_dimension =
         derived_dimension<length_base_dimension, 1,
                           mass_base_dimension, 1,
                           time_base_dimension, -2,
                           current_base_dimension, -1>::type;

      //! magnetic reluctivity dimension
      using magnetic_reluctivity_dimension =
         derived_dimension<length_base_dimension, -1,
                           mass_base_dimension, -1,
                           time_base_dimension, 2,
                           current_base_dimension, 2>::type;

      //! magnetic reluctivity derived dimension
      /*!
       *  dimension for the magnetic reluctivity derived with respect
       *  to the magnetic flux density.
       */
      using magnetic_reluctivity_derived_dimension =
         derived_dimension<length_base_dimension, -1,
                           mass_base_dimension, -2,
                           time_base_dimension, 4,
                           current_base_dimension, 3>::type;

      //! the finite difference coefficients' dimension
      using fd_coefficient_dimension =
         derived_dimension<length_base_dimension, -3,
                           mass_base_dimension, -1,
                           time_base_dimension, 2,
                           current_base_dimension, 2>::type;

      //! thermal conductivity
      using thermal_conductivity_dimension =
         derived_dimension<length_base_dimension, 1,
                           mass_base_dimension, 1,
                           time_base_dimension, -3,
                           temperature_base_dimension, -1>::type;

      namespace si {
         BOOST_UNITS_STATIC_CONSTANT(siemens_per_metre, conductivity);

         //! current density unit
         using current_density = unit<current_density_dimension, si::system>;

         //! The current density's unit
         BOOST_UNITS_STATIC_CONSTANT(ampere_per_metre_squared, current_density);

         //! magnetic vector potential unit
         using magnetic_vector_potential =
            unit<magnetic_vector_potential_dimension, si::system>;

         //! The magnetic vector potential's unit
         BOOST_UNITS_STATIC_CONSTANT(weber_per_metre,
                                     magnetic_vector_potential);

         //! magnetic reluctivity unit
         using magnetic_reluctivity =
            unit<magnetic_reluctivity_dimension, si::system>;

         //! The magnetic reluctivity's unit
         BOOST_UNITS_STATIC_CONSTANT(metre_per_henry, magnetic_reluctivity);

         //! The magnetic permeability's unit
         BOOST_UNITS_STATIC_CONSTANT(henry_per_metre, permeability);

         //! magnetic reluctivity derived unit
         using magnetic_reluctivity_derived =
            unit<magnetic_reluctivity_derived_dimension, si::system>;

         /*!
          *  Unit for the magnetic reluctivity when derived with
          *  respect to the magnetic flux density.
          */
         BOOST_UNITS_STATIC_CONSTANT(metre_per_henry_per_tesla,
                                     magnetic_reluctivity_derived);

         //! finite difference coefficients' unit
         using fd_coefficient_unit =
            unit<fd_coefficient_dimension, si::system>;

         //! The finite difference coefficient's unit
         BOOST_UNITS_STATIC_CONSTANT(ampere_per_tesla_per_metre_cubic,
                                     fd_coefficient_unit);

         //! The magnetic field intensity's unit
         BOOST_UNITS_STATIC_CONSTANT(ampere_per_metre,
                                     magnetic_field_intensity);

         using thermal_conductivity =
            unit<thermal_conductivity_dimension, si::system>;
      }

      //! Generate a in_base_units function
      /*!
       * The in_base_untis function divides a quantity by its base
       * unit.
       *
       * @param QUANTITY is a unit type.
       *
       * @param UNIT is the base unit.
       */
#define BOOST_UNITS_IN_BASE_UNITS(QUANTITY, UNIT)                       \
      template <typename Float>                                         \
      inline Float in_base_units(const quantity<QUANTITY, Float>& v)	\
      {                                                                 \
         return static_cast<Float>(v/UNIT);                             \
      }

      BOOST_UNITS_IN_BASE_UNITS(si::area, si::square_metre)
      BOOST_UNITS_IN_BASE_UNITS(si::conductivity, si::siemens_per_metre)
      BOOST_UNITS_IN_BASE_UNITS(si::current, si::ampere)
      BOOST_UNITS_IN_BASE_UNITS(si::current_density,
                                si::ampere_per_metre_squared)
      BOOST_UNITS_IN_BASE_UNITS(si::fd_coefficient_unit,
                                si::ampere_per_tesla_per_metre_cubic)
      BOOST_UNITS_IN_BASE_UNITS(si::frequency, si::hertz)
      BOOST_UNITS_IN_BASE_UNITS(si::length, si::metre)
      BOOST_UNITS_IN_BASE_UNITS(si::magnetic_field_intensity,
                                si::ampere_per_metre)
      BOOST_UNITS_IN_BASE_UNITS(si::magnetic_flux_density, si::tesla)
      BOOST_UNITS_IN_BASE_UNITS(si::magnetic_reluctivity,
                                si::metre_per_henry)
      BOOST_UNITS_IN_BASE_UNITS(si::magnetic_reluctivity_derived,
                                si::metre_per_henry_per_tesla)
      BOOST_UNITS_IN_BASE_UNITS(si::magnetic_vector_potential,
                                si::weber_per_metre)
      BOOST_UNITS_IN_BASE_UNITS(si::plane_angle, si::radians)
      BOOST_UNITS_IN_BASE_UNITS(si::time, si::second)

#ifdef BOOST_UNITS_IN_BASE_UNITS
#undef BOOST_UNITS_IN_BASE_UNITS
#endif
   }
}

#endif // BOOST_UNITS_EXTRA_HPP_FH1NGMhc_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
