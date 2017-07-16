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
 *  \brief Combine two reluctivities (useful for laminated materials)
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef COMBINED_RELUCTIVITY_HPP_tIMAHdo1_
#define COMBINED_RELUCTIVITY_HPP_tIMAHdo1_

#include <memory>
#include <tuple>

#include "config.h"

#include "scalar_reluctivity.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {

   //! Magnetic flux passing the material parallel to lamination direction
   /*!
    *  Calculate a laminated material's magnetic reluctivity and its
    *  derivation with respect to the magnetic flux density in the case
    *  the magnetic flux is passing the material parallel to the
    *  lamination direction.
    */
   struct FluxParallel {
      //! Calculate the materials magnetic reluctivity
      /*!
       *  \param reluctivity_1 is the first material's reluctivity.
       *
       *  \param reluctivity_2 is the second material's reluctivity.
       *
       *  \param fill_factor_1 is the first material's fill factor.  The
       *  second material's fill factor \f$k_2\f$ is obviously
       *  \f$k_2=1-k_1\f$.
       *
       *  \returns the laminated material's magnetic reluctivity
       */
      static MagneticReluctivity Reluctivity(
         MagneticReluctivity reluctivity_1,
         MagneticReluctivity reluctivity_2, double fill_factor_1)
      {
         return reluctivity_1 * reluctivity_2
            / Denominator(reluctivity_1, reluctivity_2, fill_factor_1);
      }

      //! Calculate the reluctivity's derivation
      /*!
       *  Calculate the material's reluctivity's derivation with respect
       *  to the magnetic flux density.
       *
       *  \param reluctivity_1 is the first material's reluctivity.
       *
       *  \param reluctivity_2 is the second material's reluctivity.
       *
       *  \param reluctivity_derived_1 is the first material's
       *  reluctivity's derivation.
       *
       *  \param reluctivity_derived_2 is the second material's
       *  reluctivity's derivation.
       *
       *  \param fill_factor_1 is the first material's fill factor.  The
       *  second material's fill factor \f$k_2\f$ is obviously
       *  \f$k_2=1-k_1\f$.
       *
       *  \returns the laminated material's reluctivity's derivation with
       *  respect to the magnetic flux density.
       */
      static MagneticReluctivityDerived ReluctivityDerived(
         MagneticReluctivity reluctivity_1,
         MagneticReluctivity reluctivity_2,
         MagneticReluctivityDerived reluctivity_derived_1,
         MagneticReluctivityDerived reluctivity_derived_2,
         double fill_factor_1)
      {
         const auto denominator = Denominator(
            reluctivity_1, reluctivity_2, fill_factor_1);

         const auto u_strich = reluctivity_derived_1 * reluctivity_2
            + reluctivity_1 * reluctivity_derived_2;

         const auto uv_strich = reluctivity_1 * reluctivity_2
            * (fill_factor_1
               * (reluctivity_derived_2 - reluctivity_derived_1)
               + reluctivity_derived_1);

         return u_strich / denominator + uv_strich / denominator / denominator;
      }

   private:
      static MagneticReluctivity Denominator(
         MagneticReluctivity reluctivity_1,
         MagneticReluctivity reluctivity_2, double fill_factor_1)
      {
         return fill_factor_1 * (reluctivity_2 - reluctivity_1) + reluctivity_1;
      }
   };

   //! Magnetic flux passing the material perpendicular to lamination direction
   /*!
    *  Calculate a laminated material's magnetic reluctivity and its
    *  derivation with respect to the magnetic flux density in the
    *  case the magnetic flux is passing the material perpendicular to
    *  the lamination direction.
    */
   struct FluxPerpendicular {
      //! Calculate the materials magnetic reluctivity
      /*!
       *  \param reluctivity_1 is the first material's reluctivity.
       *
       *  \param reluctivity_2 is the second material's reluctivity.
       *
       *  \param fill_factor_1 is the first material's fill factor.  The
       *  second material's fill factor \f$k_2\f$ is obviously
       *  \f$k_2=1-k_1\f$.
       *
       *  \returns the laminated material's magnetic reluctivity
       */
      static MagneticReluctivity Reluctivity(
         MagneticReluctivity reluctivity_1,
         MagneticReluctivity reluctivity_2, double fill_factor_1)
      {
         return fill_factor_1 * (reluctivity_1 - reluctivity_2) + reluctivity_2;
      }

      //! Calculate the reluctivity's derivation
      /*!
       *  Calculate the material's reluctivity's derivation with respect
       *  to the magnetic flux density.
       *
       *  \param reluctivity_1 is the first material's reluctivity.
       *
       *  \param reluctivity_2 is the second material's reluctivity.
       *
       *  \param reluctivity_derived_1 is the first material's
       *  reluctivity's derivation.
       *
       *  \param reluctivity_derived_2 is the second material's
       *  reluctivity's derivation.
       *
       *  \param fill_factor_1 is the first material's fill factor.  The
       *  second material's fill factor \f$k_2\f$ is obviously
       *  \f$k_2=1-k_1\f$.
       *
       *  \returns the laminated material's reluctivity's derivation with
       *  respect to the magnetic flux density.
       */
      static MagneticReluctivityDerived ReluctivityDerived(
         MagneticReluctivity reluctivity_1,
         MagneticReluctivity reluctivity_2,
         MagneticReluctivityDerived reluctivity_derived_1,
         MagneticReluctivityDerived reluctivity_derived_2,
         double fill_factor_1)
      {
         return fill_factor_1 * (reluctivity_derived_1 - reluctivity_derived_2)
            + reluctivity_derived_2;
      }
   };

   //! A laminated material's reluctivity
   /*!
    *  \tparam Orientation describes whether the flux is passing the
    *  material parallel or perpendicular to the lamination direction.
    */
   template<typename Orientation>
   struct CombinedReluctivity:
      ScalarReluctivity,
      detail::MakeShared<CombinedReluctivity<Orientation> >
   {
   private:
      using ReluctivityPointer = std::shared_ptr<const ScalarReluctivity>;

   public:
      //! Create a combined reluctivity
      /*!
       *  \param reluctivity_1 is the first layer's reluctivity
       *
       *  \param reluctivity_2 is the second layer's reluctivity
       *
       *  \param fill_factor \f$k_{\mathrm{f}}\f$ is the fill factor
       *  related to reluctivity_1:
       *  \f[k_{\mathrm{f}}=\frac{d_1}{d_1+d_2}\f]
       */
      CombinedReluctivity(
         ReluctivityPointer reluctivity_1,
         ReluctivityPointer reluctivity_2, double fill_factor)
         : reluctivity_1_(reluctivity_1), reluctivity_2_(reluctivity_2),
         fill_factor_(fill_factor)
      {
      }

      ReluctivityReluctivityDerivedPair operator ()(
         MagneticFluxDensityReal b) const override;

      bool ResultCacheRecommended() const override
      {
         return reluctivity_1_->ResultCacheRecommended()
            or reluctivity_2_->ResultCacheRecommended();
      }

      //! Destruction
      virtual ~CombinedReluctivity()
      {
      }

   private:
      ReluctivityPointer reluctivity_1_;
      ReluctivityPointer reluctivity_2_;
      double fill_factor_;
   };

   template<typename Orientation>
   ReluctivityReluctivityDerivedPair
   CombinedReluctivity<Orientation>::operator ()(
      MagneticFluxDensityReal b) const
   {
      using std::get;

      ReluctivityReluctivityDerivedPair
         pair_1 = (*reluctivity_1_)(b), pair_2 = (*reluctivity_2_)(b);

      return std::make_tuple(
         Orientation::Reluctivity(
            get<0>(pair_1), get<0>(pair_2), fill_factor_),
         Orientation::ReluctivityDerived(
            get<0>(pair_1), get<0>(pair_2),
            get<1>(pair_1), get<1>(pair_2), fill_factor_));
   }

#ifndef INSTANTIATE_
   extern template struct CombinedReluctivity<FluxParallel>;
   extern template struct CombinedReluctivity<FluxPerpendicular>;
#endif

}

#endif // COMBINED_RELUCTIVITY_HPP_tIMAHdo1_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
