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
 *  \brief edge properties
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef EDGE_PROPERTIES_HPP_o6p09D7f_
#define EDGE_PROPERTIES_HPP_o6p09D7f_

#include <cmath>
#include <iostream>

#include "config.h"

#include <boost/format.hpp>
#include <boost/units/cmath.hpp>

#include "analysis.hpp"
#include "coordinate_systems.hpp"
#include "exception.hpp"
#include "grid_graph-fwd.hpp"
#include "scalar_reluctivity.hpp"
#include "utility.hpp"

#ifdef HAVE_LIBLOKI
#include <loki/SmallObj.h>
#endif

namespace PACKAGE_NS {
   namespace fi
   {
      struct ReluctivityAccess {
         virtual void InvalidateCache () const =0;
         virtual ReluctivityReluctivityDerivedPair operator() (
            MagneticFluxDensityReal) const =0;
      };

      struct CachedReluctivity : ReluctivityAccess
      {
      public:
         CachedReluctivity(ConstReluctivityPointer reluctivity)
            : reluctivity_(reluctivity), cache_valid_(false)
         {
         }

         void InvalidateCache() const override final {
            cache_valid_ = false;
         }

         ReluctivityReluctivityDerivedPair operator() (
            MagneticFluxDensityReal flux_density) const override final
         {
            if (not cache_valid_) {
               cached_reluctivity_ = (*reluctivity_)(flux_density);
               cache_valid_ = true;
            }

            return cached_reluctivity_;
         }

      private:
         ConstReluctivityPointer reluctivity_;
         mutable bool cache_valid_;
         mutable ReluctivityReluctivityDerivedPair cached_reluctivity_;
      };

      struct DirectReluctivity : ReluctivityAccess
      {
         DirectReluctivity(ConstReluctivityPointer reluctivity)
            : reluctivity_(reluctivity)
         {
         }

         void InvalidateCache() const override final
         {
         }

         ReluctivityReluctivityDerivedPair operator() (
            MagneticFluxDensityReal flux_density) const override final
         {
            return (*reluctivity_)(flux_density);
         }

      private:
         ConstReluctivityPointer reluctivity_;
      };

      //! Grid edge properties
      /*!
       * @tparam Analysis is the analysis type
       */
      template <typename Analysis>
      struct GridEdge
      {
      private:
         template <typename T>
         using PassType = ::PACKAGE_NS::detail::PassType<T>;
         using CoefficientQuantity =
            boost::units::quantity<boost::units::si::magnetic_reluctivity,
                                   typename Analysis::DataType>;
         using CoefficientPair =
            std::tuple<CoefficientQuantity, CoefficientQuantity>;
         using MagneticFieldIntensity =
            typename Analysis::MagneticFieldIntensity;
         using MagneticFluxDensity = typename Analysis::MagneticFluxDensity;
         using MagneticVoltage = typename Analysis::MagneticVoltage;
         using Point = CartesianSystem::Point;
         using ReluctivityAccessPointer =
            std::unique_ptr<ReluctivityAccess>;

         static ReluctivityAccessPointer MakeReluctivityHolder(
            bool cached, ConstReluctivityPointer reluctivity)
         {
            if (cached)
               return std::make_unique<CachedReluctivity>(reluctivity);

            return std::make_unique<DirectReluctivity>(reluctivity);
         }

      public:
         //! Initialization
         /*!
          * @param edge_length the edge'es length
          *
          * @param length_coefficient length of the secondary grid's
          * crossing edge devided by the edge'es length.
          *
          * @param reluctivity reluctivity model for the direction
          * perpendicular to the edge.
          */
         GridEdge(Length edge_length, double length_coefficient,
                  ConstReluctivityPointer reluctivity)
            : edge_length_(edge_length),
              length_coefficient_(length_coefficient),
              reluctivity_holder_(
                 MakeReluctivityHolder(
                    reluctivity->ResultCacheRecommended(), reluctivity)),
              flux_density_(MagneticFluxDensity::from_value(.0))
         {
            this->CheckValid(length_coefficient_);
         }

         ReluctivityReluctivityDerivedPair EvalReluctivity() const
         {
            return (*reluctivity_holder_)(
               boost::units::abs(this->flux_density_));
         }

         //! Compute the reluctivity at this edge
         MagneticReluctivity reluctivity() const
         {
            return std::get<0>(this->EvalReluctivity());
         }

         //! Compute the system matrix coefficient for this edge
         CoefficientQuantity SystemMatrixCoefficient() const
         {
            return reluctivity() * length_coefficient_;
         }

         //! Compute the system and jacobi matrix coefficient
         CoefficientPair SystemAndJacobiMatrixCoefficient() const
         {
            using std::get;

            const auto reluctivity_reluctivity_derived
               = this->EvalReluctivity();

            const CoefficientQuantity std_coeff =
               length_coefficient_ * get<0>(reluctivity_reluctivity_derived);

            const CoefficientQuantity jac_coeff = std_coeff
               + length_coefficient_ * abs(flux_density_)
               * std::get<1>(reluctivity_reluctivity_derived);

            CheckValid(std_coeff);
            CheckValid(jac_coeff);

            return std::make_tuple(std_coeff, jac_coeff);
         }

         Length length() const {
            return edge_length_;
         }

         const double& length_coefficient() const {
            return length_coefficient_;
         }

         Length secondary_length() const {
            return length_coefficient_*edge_length_;
         }

         MagneticFluxDensity flux_density() const {
            return flux_density_;
         }

         void set_flux_density(PassType<MagneticFluxDensity> new_flux_density)
         {
            flux_density_ = new_flux_density;
            reluctivity_holder_->InvalidateCache();
         }

         MagneticFieldIntensity magnetic_field_intensity() const {
            return this->reluctivity() * this->flux_density_;
         }

         MagneticVoltage magnetic_voltage() const {
            return this->magnetic_field_intensity()
               * secondary_length();
         }

      private:
         void CheckValid(const double& v) const
         {
            if (not std::isfinite(v))
               throw RuntimeError("Invalid length coefficient!");
         }

         void CheckValid(const CoefficientQuantity& c) const
         {
#ifndef NDEBUG
            if (not std::isfinite(in_base_units(abs(c))))
               throw RuntimeError("Invalid matrix coefficient!");
#endif
         }

      private:
      	 Length edge_length_;
         double length_coefficient_;
         ReluctivityAccessPointer reluctivity_holder_;
         MagneticFluxDensity flux_density_;
      };

      template <typename Analysis>
      struct StraightGridEdge : GridEdge<Analysis>
      {
      private:
         using Base = GridEdge<Analysis>;

      public:
         StraightGridEdge(Length edge_length, double length_coefficient,
                          ConstReluctivityPointer reluctivity)
            : Base(edge_length, length_coefficient, reluctivity)
         {
         }
      };

      template <typename Analysis>
      struct CircularGridEdge : GridEdge<Analysis>
      {
      private:
         using Base = GridEdge<Analysis>;

      public:
         CircularGridEdge(Length edge_length, double length_coefficient,
                          Length radius,
                          ConstReluctivityPointer reluctivity)
            : Base(edge_length, length_coefficient, reluctivity),
              radius_(radius)
         {
         }

      private:
         Length radius_;
      };

#ifndef INSTANTIATE_
      extern template struct StraightGridEdge<MagnetoHarmonicAnalysis<> >;
      extern template struct StraightGridEdge<TimeTransientAnalysis<> >;
      extern template struct CircularGridEdge<MagnetoHarmonicAnalysis<> >;
      extern template struct CircularGridEdge<TimeTransientAnalysis<> >;
#endif
   }
}

#endif // EDGE_PROPERTIES_HPP_o6p09D7f_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
