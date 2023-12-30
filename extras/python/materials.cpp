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

#include <boost/parameter/python.hpp>

#include "config.h"

#include "boost_python_extra.hpp"

#include "coordinate_systems.hpp"
#include "material.hpp"
#include "reluctivity_models.hpp"

namespace PACKAGE_NS {
   namespace python {

      // *** ArctanReluctivity ***

      struct ArctanReluctivityWrap : ArctanReluctivity
      {
      private:
         using Base = ArctanReluctivity;

      public:
         ArctanReluctivityWrap(
            double initial_mu_r,
            MagneticPolarizationReal saturation_polarization) :
            Base(initial_mu_r, saturation_polarization)
         {
         }

         ArctanReluctivityWrap(
            bp::object h_values, bp::object bi_values) :
            Base(
               bp::stl_input_iterator<MagneticFieldIntensityReal>(h_values),
               bp::stl_input_iterator<MagneticFieldIntensityReal>(),
               bp::stl_input_iterator<MagneticPolarizationReal>(bi_values))
         {
         }
      };


      // *** SqrtReluctivity ***

      struct SqrtReluctivityWrap : SqrtReluctivity
      {
      private:
         using Base = SqrtReluctivity;

      public:
         SqrtReluctivityWrap(
            double initial_mu_r,
            MagneticPolarizationReal saturation_polarization,
            double a) : Base(
               initial_mu_r,
               saturation_polarization,
               a)
         {
         }

         SqrtReluctivityWrap(bp::object h_values, bp::object bi_values)
            : Base(
               bp::stl_input_iterator<MagneticFieldIntensityReal>(h_values),
               bp::stl_input_iterator<MagneticFieldIntensityReal>(),
               bp::stl_input_iterator<MagneticPolarizationReal>(bi_values))
         {
         }
      };


      // *** SplineReluctivity ***

      struct SplineReluctivityWrap : SplineReluctivity
      {
      private:
         using Base = SplineReluctivity;

      public:
         SplineReluctivityWrap(bp::object h_values, bp::object bi_values)
            : Base(
               bp::stl_input_iterator<MagneticFieldIntensityReal>(h_values),
               bp::stl_input_iterator<MagneticFieldIntensityReal>(),
               bp::stl_input_iterator<MagneticPolarizationReal>(bi_values))
         {
         }
      };


      // *** CombinedReluctivity ***

      template <typename Orientation>
      void DefCombinedReluctivity(const char* name)
      {
         namespace bp = boost::python;

         BPClassSharedPtr<CombinedReluctivity<Orientation> >(
            name, bp::no_init)
            .def(bp::init<ConstReluctivityPointer,
                 ConstReluctivityPointer, double>())
            .def("__call__", &ArctanReluctivityWrap::operator());
      }

      template <typename ReluctivityModel>
      void AddImplicitPointerConversion(ReluctivityModel*)
      {
         boost::python::implicitly_convertible<
            std::shared_ptr<ReluctivityModel>,
            PACKAGE_NS::ConstReluctivityPointer>();
      }
   }
}

BOOST_PYTHON_MODULE(MODULE_NAME)
{
   using namespace PACKAGE_NS;
   using namespace PACKAGE_NS::python;
   namespace bp = boost::python;
   namespace mpl = boost::mpl;
   namespace bpp = boost::parameter::python;

   using ComboReluctivityFluxPar = CombinedReluctivity<FluxParallel>;
   using ComboReluctivityFluxPerp = CombinedReluctivity<FluxPerpendicular>;

   // *** ArctanReluctivity ***
   BPClassSharedPtr<ArctanReluctivityWrap>(
      "ArctanReluctivity", bp::no_init)
      .def(bp::init<bp::object, bp::object>())
      .def(bp::init<double, MagneticPolarizationReal>())
      .def("__call__", &ArctanReluctivityWrap::operator())
      .add_property(
         "initial_mu_r", &ArctanReluctivityWrap::initial_mu_r)
      .add_property("saturation_polarization",
                    &ArctanReluctivityWrap::saturation_polarization);

   // *** SqrtReluctivity ***
   BPClassSharedPtr<SqrtReluctivityWrap>("SqrtReluctivity", bp::no_init)
      .def(bp::init<bp::object, bp::object>())
      .def(bp::init<double, MagneticPolarizationReal, double>())
      .def("__call__", &SqrtReluctivityWrap::operator())
      .add_property("initial_mu_r", &SqrtReluctivityWrap::initial_mu_r)
      .add_property("saturation_polarization",
                    &SqrtReluctivityWrap::saturation_polarization)
      .add_property("a", &SqrtReluctivityWrap::a);

   // *** SplineReluctivity ***
   BPClassSharedPtr<SplineReluctivityWrap>("SplineReluctivity", bp::no_init)
      .def(bp::init<bp::object, bp::object>())
      .def("__call__", &SplineReluctivityWrap::operator());

   // *** ConstantReluctivity ***
   BPClassSharedPtr<ConstantReluctivity>("ConstantReluctivity", bp::no_init)
      .def(bp::init<double>())
      .def(bp::init<Permeability>())
      .def(bp::init<MagneticReluctivity>())
      .def("__call__", &ConstantReluctivity::operator());

   // *** CombinedReluctivity ***
   DefCombinedReluctivity<FluxParallel>("ComboReluctivityFluxPar");
   DefCombinedReluctivity<FluxPerpendicular>("ComboReluctivityFluxPerp");

   // *** Material ***
   BPClassSharedPtr<Material>("Material", bp::no_init)
      .def(bpp::init< mpl::vector<
           tag::anisotropic_reluctivity*(ReluctivityTuple),
           tag::isotropic_reluctivity*(ConstReluctivityPointer),
           tag::electrical_conductivity*(ElectricalConductivity),
           tag::thermal_conductivity*(ThermalConductivity),
           tag::permittivity*(Permittivity)> >() )
      .def("EvaluateConductivity", &Material::EvaluateConductivity)
      .def("EvaluateReluctivity", &Material::EvaluateReluctivity);

   bp::enum_<Axis>("Axis")
      .value("first", Axis::kFirst)
      .value("second", Axis::kSecond);

   AddImplicitPointerConversion((ArctanReluctivityWrap*)nullptr);
   AddImplicitPointerConversion((ComboReluctivityFluxPerp*)nullptr);
   AddImplicitPointerConversion((ComboReluctivityFluxPar*)nullptr);
   AddImplicitPointerConversion((ConstantReluctivity*)nullptr);
   AddImplicitPointerConversion((SplineReluctivityWrap*)nullptr);
   AddImplicitPointerConversion((SqrtReluctivityWrap*)nullptr);

   bp::implicitly_convertible<
      std::shared_ptr<Material>, ConstMaterialPointer>();
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
