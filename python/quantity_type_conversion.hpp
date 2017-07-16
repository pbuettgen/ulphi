// -*- coding: utf-8 -*-
/*
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and library
 * for computing electromagnetic fields.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*! \file
 *
 *  \brief Conversion between boost::units::quantity and Python Pint quantities
 *  \author Philipp Büttgenbach
 *  \date $Date: 2016-02-14 23:53:41 +0100 (So, 14. Feb 2016) $
 *  \copyright MPL 2.0
 *
 * Conversion between
 *   <a href="http://localhost/cgi-bin/dwww/usr/share/doc/libboost1.49-doc/HTML/doc/html/boost_units.html">boost::units::quantity</a> and <a href="http://localhost/cgi-bin/dwww/usr/share/doc/python-pint-doc/html/index.html?type=html">Python Pint</a> quantities.
 */

#ifndef QUANTITY_CONVERSION_HPP_ZZQgIzgw_
#define QUANTITY_CONVERSION_HPP_ZZQgIzgw_

#include <cstring>

#include <boost/units/io.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/python.hpp>

#include "config.h"
#include "boost_units_extra.hpp"

namespace PACKAGE_NS
{
   namespace python
   {
      //! Convert boost quantities to python pint quantities
      /*!
       * \tparam Unit is the unit type
       * \tparam Y is the numeric type
       */
      template<class Unit, class Y = double>
      struct QuantityToPythonConversion
      {
         static PyObject* convert(
            const boost::units::quantity<Unit, Y> & q)
         {
            namespace bp = boost::python;

            static bp::object qty_class = ImportQty();

            bp::object pint_qty
               = qty_class(q.value(), symbol_string(Unit()));

            return bp::incref(pint_qty.ptr());
         }

      private:
         static boost::python::object ImportQty() {
            namespace bp = boost::python;

            const bp::object ulphi_quantity_module
               = bp::import("ulphi._quantity");

            return ulphi_quantity_module.attr("Qty");
         }
      };


      namespace detail
      {
#define PYTHON_PINT_DIMENSION(UNIT, UNITSTR)            \
         constexpr const char* PythonPintDimension (    \
            boost::units::UNIT* ) {                     \
            return UNITSTR;                             \
         }

         PYTHON_PINT_DIMENSION(si::conductivity, "[current] ** 2 * [time] ** 3 / [length] ** 3 / [mass]")
         PYTHON_PINT_DIMENSION(si::current, "[current]")
         PYTHON_PINT_DIMENSION(si::current_density, "[current] / [length] ** 2")
         PYTHON_PINT_DIMENSION(si::frequency, "1 / [time]")
         PYTHON_PINT_DIMENSION(si::length, "[length]")
         PYTHON_PINT_DIMENSION(si::magnetic_field_intensity, "[current] / [length]")
         PYTHON_PINT_DIMENSION(si::magnetic_flux_density, "[mass] / [current] / [time] ** 2")
         PYTHON_PINT_DIMENSION(si::magnetic_reluctivity, "[current] ** 2 * [time] ** 2 / [length] / [mass]")
         PYTHON_PINT_DIMENSION(si::magnetic_reluctivity_derived, "[current] ** 3 * [time] ** 4 / [length] / [mass] ** 2")
         PYTHON_PINT_DIMENSION(si::magnetic_vector_potential, "[length] * [mass] / [current] / [time] ** 2")
         PYTHON_PINT_DIMENSION(si::permeability, "[length] * [mass] / [current] ** 2 / [time] ** 2")
         PYTHON_PINT_DIMENSION(si::plane_angle, "dimensionless")
         PYTHON_PINT_DIMENSION(si::time, "[time]")

#undef PYTHON_PINT_DIMENSION
      }

      //! Check whether a python object can be converted to a quantity
      /*!
       * A python object is convertible if
       *
       * - it has a magnitude property
       * - the magnitude's type matches
       * - it has a unit property
       * - it has a dimensionality property
       * - the dimensionality matches
       *
       * \tparam Unit is the unit type
       * \tparam Y is the numeric type
       */
      template <typename Unit, typename Y = double>
      void * QuantityIsConvertible ( PyObject* o )
      {
         namespace bp = boost::python;

         bp::expect_non_null(o);

         if (not PyObject_HasAttrString(o, "magnitude"))
            return nullptr;

         // Check magnitude
         bp::object magnitude (
            bp::handle<> (
               PyObject_GetAttrString ( o, "magnitude" ) ) );

         bp::extract<Y> magnitude_extract ( magnitude );

         if ( not magnitude_extract.check() ) return nullptr;

         // Check units
         if ( not PyObject_HasAttrString ( o, "units" ) )
            return nullptr;

         // Check dimensionality
         if (not PyObject_HasAttrString(o, "dimensionality"))
            return nullptr;

         bp::object dimensionality (
            bp::handle<> (
               PyObject_GetAttrString ( o, "dimensionality" ) ) );

         const std::string dimensionality_string
            = bp::extract<std::string>(bp::str(dimensionality));

         if (dimensionality_string != detail::PythonPintDimension(
                static_cast<Unit*>(nullptr)))
            return nullptr;

         return o;
      }

      namespace detail {
         template<typename Unit, typename Y = double>
         boost::units::quantity<Unit, Y> ConvertFromPythonPint(
            PyObject* o)
         {
            namespace bp = boost::python;
            namespace bu = boost::units;

            using Quantity = bu::quantity<Unit, Y>;

            bp::object pint_quantity_with_unit_adjusted(
               bp::handle<>(
                  PyObject_CallMethod(
                     o, "to", "s", bu::symbol_string(Unit()).c_str() )));

            auto pint_quantity_magnitude
               = pint_quantity_with_unit_adjusted.attr ( "magnitude" );
            Y magnitude = bp::extract<Y> (pint_quantity_magnitude);

            return Quantity::from_value(magnitude);
         }
      }

      //! Convert a python object to a boost::units::quantity object
      /*!
       * \tparam Unit is the unit type
       * \tparam Y is the numeric type
       */
      template <typename Unit, typename Y = double>
      void QuantityConstructFromPint (
         PyObject* o,
         boost::python::converter::rvalue_from_python_stage1_data* data )
      {
         namespace bp = boost::python;
         namespace bu = boost::units;

         using Quantity = bu::quantity<Unit, Y>;

         // Grab pointer to memory into which to construct the new Quantity
         void * storage
            = reinterpret_cast<
            bp::converter::rvalue_from_python_storage<Quantity>* > (
               data )->storage.bytes;

         const Quantity result
            = detail::ConvertFromPythonPint<Unit, Y>(o);

         std::memcpy(storage, &result, sizeof(Quantity));

         //new ( storage ) Quantity ( Quantity::from_value ( magnitude ) );

         data->convertible = storage;
      }

      //! Register to and from python conversion for boost::units::quantity
      /*!
       * \tparam Unit is the unit type
       * \tparam Y is the numeric type
       */
      template <typename Unit, typename Y = double>
      void QuantityRegisterConversion(
         boost::units::quantity<Unit, Y>* = nullptr)
      {
         namespace bp = boost::python;

         using Quantity = boost::units::quantity<Unit, Y>;

         // register the to-python converter
         bp::to_python_converter<
            Quantity, QuantityToPythonConversion<Unit, Y> >();

         // register the from-python converter
         bp::converter::registry::push_back (
            &QuantityIsConvertible<Unit, Y>,
            &QuantityConstructFromPint<Unit, Y>,
            bp::type_id<Quantity>() );
      }
   }
}

#endif

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
