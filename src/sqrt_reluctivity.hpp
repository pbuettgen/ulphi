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
 *  \brief Use a sqrt-function to model magnetic reluctivity behaviour
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef SQRT_RELUCTIVITY_HPP_yyyoft3A_
#define SQRT_RELUCTIVITY_HPP_yyyoft3A_

#include <map>

#include "config.h"

#include <boost/units/systems/si/codata_constants.hpp>

#include "eigen_extra.hpp"
#include "lm_functor.hpp"
#include "scalar_reluctivity.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Fuctor used by the Levenberg-Marquardt fitting procedure
      struct LMFunctorSqrt: LMFunctor<3> {
      private:
         using Base = LMFunctor<3>;

      public:
         //! Initialization
         /*!
          * @tparam HIterator is a type for iterating over magnetic
          * field intensity values.
          *
          * @tparam BiIterator is a type for iterating over magnetic
          * polarization values.
          *
          * @param h_begin points to the start of a sequence with
          * magnetic field intensity values
          *
          * @param h_end marks the end of the field intensity value's sequence
          *
          * @param bi_begin points to the start of a sequence with
          * magnetic polarization values.
          */
         template<typename HIterator, typename BiIterator>
         LMFunctorSqrt(HIterator h_begin, const HIterator& h_end,
                       BiIterator bi_begin) :
            Base(h_begin, h_end, bi_begin)
         {
         }

         //! Compute delta for current Levenberg-Marquardt step
         /*!
          * @param [in]  x is the current parameter vector
          * @param [out] fvec is for storing the computed deltas
          * @return status/error code
          */
         int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
         {
            using boost::units::si::constants::codata::mu_0;
            const double mu_0_value = mu_0.value().value();

            // x(0) -> intitial mu_r
            // x(1) -> saturation polarization Js
            // x(2) -> a

            const auto a2 = 1. - x(2);
            const Eigen::ArrayXd Ha = mu_0_value * bi_h_values_.col(0)
               * (x(0) - 1.) / x(1);

            fvec = x(1) * (Ha + 1. - sqrt(Pow<2>(Ha + 1.) - 4. * Ha * a2))
               / (2. * a2) - bi_h_values_.col(1);

            return 0;
         }

         //! Compute the jacobi matrix
         /*!
          * @param [in]  x is the actual parameter vector
          * @param [out] fjac is the jacobi matrix
          * @return status/error code
          */
         int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const
         {
            using boost::units::si::constants::codata::mu_0;

            const auto& mu_r = x(0);
            const auto& Js = x(1);
            const auto& a = x(2);
            const auto& H = bi_h_values_.col(0);
            const double mu0 = mu_0.value().value();

            // --------------------------------------------------------

            const auto a2 = 1 - a;
            const auto a3 = 2 * a2;
            const auto h01 = H * mu0 / Js;
            const auto h02 = h01 * (mu_r - 1) + 1;

            fjac.col(0) = Js
               * (h01
                  - h01 * (-2 * a2 + h02)
                  / sqrt(
                     -4 * h01 * (-a + 1) * (mu_r - 1)
                     + Pow<2>(h02))) / a3;

            // --------------------------------------------------------

            const auto h11 = H * mu0 * (mu_r - 1);
            const auto h12 = sqrt(
               (-4 * h11 * Js * a2 + Pow<2>(h11 + Js)) / Pow<2>(Js));

            fjac.col(1) = (h11 * (h11 - Js * h12 + (-a3 + 1) * Js)
                           + Js * h12 * (h11 - Js * h12 + Js)) / (Pow<2>(Js) * h12 * a3);

            // --------------------------------------------------------
            const auto h21 = H * mu0 * (mu_r - 1);
            const auto h22 = h21 / Js;
            const auto h23 = sqrt(-4 * h22 * a2 + Pow<2>(h22 + 1));

            fjac.col(2) = -2 * h21 / (a3 * h23)
               + 2 * Js * (h22 - h23 + 1) / Pow<2>(a3);

            return 0;
         }
      };
   }

   //! Use a sqrt-function to model magnetic reluctivity behaviour
   /*!
    * \f{align*}{
    * B(H)&=\mu_{0}H + J_{\mathrm{s}}\frac{H_{\mathrm{a}}+1-sqrt{(H_{\mathrm{a}}+1)^2-4H_{\mathrm{a}}(1-a)}}{2(1-a)}
    * \intertext{with}
    * H_{\mathrm{a}}&=\mu_0H\frac{mu_{\mathrm{r}}-1}{J_{\mathrm{s}}} \\
    * \mu_{\mathrm{r}} &= \text{initial relative permeability} \\
    * J_{\mathrm{s}} &= \text{saturation polarization} \\
    * a &= \text{knee adjust parameter}
    * \f}
    */
   struct SqrtReluctivity: ScalarReluctivity {
   private:
      template<typename T>
      using AddConstAndRef = detail::AddConstAndRef<T>;

   public:
      //! Initialization
      /*!
       * @param initial_mu_r is the relative permeability at H=0 A/m
       * @param saturation_polarization is the saturation polarization
       * @param a is the knee adjust parameter
       */
      SqrtReluctivity(double initial_mu_r,
                      MagneticPolarizationReal saturation_polarization, double a) :
         initial_mu_r_(initial_mu_r), saturation_polarization_(
            saturation_polarization), a_(a)
      {
      }

      SqrtReluctivity(const SqrtReluctivity& source) :
         initial_mu_r_(source.initial_mu_r_), saturation_polarization_(
            source.saturation_polarization_), a_(source.a_)
      {
      }

      SqrtReluctivity& operator =(const SqrtReluctivity& source) {
         initial_mu_r_ = source.initial_mu_r_;
         saturation_polarization_ = source.saturation_polarization_;
         a_ = source.a_;

         return *this;
      }

      //! Determine parameters by fitting the curve to given data points
      /*!
       * @tparam HIterator is a type for iterating over magnetic field
       * intensity values
       *
       * @tparam BiIterator is a type for iterating over magnetic
       * polarization values
       *
       * @param h_begin points to the start of a sequence with
       * magnetic field intendity values
       *
       * @param h_end points just behind this sequence
       *
       * @param bi_begin points to the start of a sequence with
       * magnetic polarization values
       */
      template<typename HIterator, typename BiIterator>
      SqrtReluctivity(HIterator, AddConstAndRef<HIterator>,
                      BiIterator);

   protected:
      //! Initialization
      template <typename HIterator, typename BiIterator>
      void Init(HIterator, AddConstAndRef<HIterator>,
                BiIterator);

   public:
      virtual ReluctivityReluctivityDerivedPair operator()(
         MagneticFluxDensityReal b) const override;

      bool ResultCacheRecommended() const override final;

   private:
      ReluctivityReluctivityDerivedPair CalculateWithNewton(
         MagneticFluxDensityReal b) const;

      ReluctivityReluctivityDerivedPair CalculateDirectly(
         MagneticFluxDensityReal b) const;

   public:
      template <typename... Args>
      static inline std::shared_ptr<SqrtReluctivity> New(Args... arguments) {
         return std::make_shared<SqrtReluctivity>(arguments...);
      }

      //! Retrieve the initial mu_r coefficient
      double initial_mu_r() const
      {
         return initial_mu_r_;
      }

      //! Retrieve the saturation polarization
      MagneticPolarizationReal saturation_polarization() const
      {
         return saturation_polarization_;
      }

      //! Retrieve parameter a
      double a() const
      {
         return a_;
      }

   private:
      using DoubleTuple = std::tuple<double, double>;

      //! Helper function to solve the B-H-curve-equation for nu
      /*!
       * Since the arctan-magnetization-curve can not be solved
       * analyticaly, it is solved numericaly using the newton
       * procedure.
       *
       * @param b is the magnetic flux densitiy under consideration
       * @param nu is the reluctivity in the current newton step
       */
      DoubleTuple NewtonIteration(MagneticFluxDensityReal, double) const;

   private:
      double initial_mu_r_;
      MagneticPolarizationReal saturation_polarization_;
      double a_;
      mutable std::map<MagneticFluxDensityReal,
                       MagneticReluctivity> value_cache_;
   };

   template <typename HIterator, typename BiIterator>
   SqrtReluctivity::SqrtReluctivity(
      HIterator h_begin, AddConstAndRef<HIterator> h_end,
      BiIterator bi_begin)
   {
      this->Init(h_begin, h_end, bi_begin);
   }

   template<typename HIterator, typename BiIterator>
   void SqrtReluctivity::Init(
      HIterator h_begin, AddConstAndRef<HIterator> h_end,
      BiIterator bi_begin)
   {
      using boost::units::si::tesla;

      Eigen::VectorXd x(3);
      // First try
      x << 5e3, 1.8, .1;
      detail::LMFunctorSqrt lmfunc(h_begin, h_end, bi_begin);
      Eigen::LevenbergMarquardt<detail::LMFunctorSqrt> lm(lmfunc);
      const auto info = lm.minimize(x);
      Eigen::CheckLevenbergMarquardtStatus(info);

      initial_mu_r_ = x(0);
      saturation_polarization_ = x(1) * tesla;
      a_ = x(2);
   }
}

#endif // SQRT_RELUCTIVITY_HPP_yyyoft3A_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
