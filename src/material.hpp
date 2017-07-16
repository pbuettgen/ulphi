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
 *  \brief Describe material properties
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef MATERIAL_HPP_v1YCe06b_
#define MATERIAL_HPP_v1YCe06b_

#include <functional>
#include <memory>
#include <string>

#include <boost/parameter.hpp>

#include "config.h"

#include "constant_reluctivity.hpp"
#include "coordinate_systems.hpp"
#include "material-fwd.hpp"
#include "scalar_reluctivity.hpp"
#include "types.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {

   BOOST_PARAMETER_NAME(anisotropic_reluctivity)
   BOOST_PARAMETER_NAME(electrical_conductivity)
   BOOST_PARAMETER_NAME(isotropic_reluctivity)
   BOOST_PARAMETER_NAME(permittivity)
   BOOST_PARAMETER_NAME(thermal_conductivity)

   using ReluctivityTuple=std::tuple<
      ConstReluctivityPointer, ConstReluctivityPointer>;

   namespace detail {
      inline ReluctivityTuple MakeFromIsotropic(
         ConstReluctivityPointer reluctivity) {
         return std::make_tuple(reluctivity, reluctivity);
      }

      struct FromIsotropic {
         FromIsotropic(ConstReluctivityPointer reluctivity)
            : reluctivity_(reluctivity)
         {}

         ReluctivityTuple operator() () const {
            return MakeFromIsotropic(reluctivity_);
         }

      private:
         ConstReluctivityPointer reluctivity_;
      };
   }

   struct MaterialBase
   {
      template <class ArgumentPack>
      MaterialBase(const ArgumentPack& args)
         : reluctivity_(
            args[_anisotropic_reluctivity
                 | detail::MakeFromIsotropic(
                    args[_isotropic_reluctivity
                         | ConstantReluctivity::New(1.) ] ) ] ),
           electrical_conductivity_(
              args[_electrical_conductivity
                   | ElectricalConductivity::from_value(.0) ] ),
           permittivity_(
              args[_permittivity
                   | boost::units::si::constants::codata::epsilon_0 ] ),
           thermal_conductivity_(
              args[_thermal_conductivity
                   | ThermalConductivity::from_value(0.) ] )
      {}

      //! Compute the material's reluctivity
      /*!
       * @param b is the magnetic flux density corresponding to the
       * requested reluctivity value.
       *
       * @return the material's reluctivity
       */
      ReluctivityReluctivityDerivedPair EvaluateReluctivity(
         MagneticFluxDensityReal b, Axis axis) const
      {
         return Axis::kFirst == axis ?
            std::get<0>(reluctivity_)->operator()(b) :
            std::get<1>(reluctivity_)->operator()(b);
      }

      //! Get a reference to the reluctivity model in use
      /*!
       * @param axis is the axis for which the reluctivity shall be
       * retrieved.
       *
       * @return a shared_ptr to the reluctivity model used.
       */
      ConstReluctivityPointer GetReluctivityModel(Axis axis) const
      {
         return Axis::kFirst == axis ?
            std::get<0>(reluctivity_) : std::get<1>(reluctivity_);
      }

      //! Compute the material's conductivity
      /*!
       * @return the material's conductivity
       */
      ElectricalConductivity EvaluateConductivity() const
      {
         return electrical_conductivity_;
      }

   private:
      ReluctivityTuple reluctivity_;
      ElectricalConductivity electrical_conductivity_;
      Permittivity permittivity_;
      ThermalConductivity thermal_conductivity_;
   };

   //! Material properties
   struct Material : MaterialBase, detail::MakeShared<Material>
   {
      //! Initialize an isotropic material with no conductivity
      /*!
       * @param reluctivity is the material's reluctivity
       */
      BOOST_PARAMETER_CONSTRUCTOR(
         Material, (MaterialBase), tag,
         (optional
          (anisotropic_reluctivity, *)
          (isotropic_reluctivity, *)
          (electrical_conductivity, *)
          (thermal_conductivity, *)
          (permittivity, *)))

      Material(const Material&) = delete;
   };
}

#endif // MATERIAL_HPP_v1YCe06b_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
