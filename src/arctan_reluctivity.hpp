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
 *  \brief Describe magnetic reluctivity using an arctan-model
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef ARCTAN_RELUCTIVITY_HPP_DvG5eh2X_
#define ARCTAN_RELUCTIVITY_HPP_DvG5eh2X_

#include <map>

#include "config.h"

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/units/systems/si/codata_constants.hpp>

#include "eigen_extra.hpp"
#include "lm_functor.hpp"
#include "scalar_reluctivity.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Functor required by the Levenberg-Marquardt procedure
      struct LMFunctorArctan: LMFunctor<2> {
         using Base = LMFunctor<2>;

         //! pi
         static double Pi() {
            return boost::math::double_constants::pi;
         }

         //! Initialization
         /*!
          * @tparam HIterator is a type for iterating over magnetic
          * field intensity values
          *
          * @tparam BiIterator is a type for iterating over magnetic
          * polarization values
          *
          * @param h_begin points to the start of magnetic field
          * intensity values.
          *
          * @param h_end points just beyond the end of these values
          *
          * @param bi_begin points to the start of the corresponding
          * polarization values.
          */
         template<typename HIterator, typename BiIterator>
         LMFunctorArctan(HIterator h_begin, const HIterator& h_end,
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
            using Eigen::atan;
            using boost::units::si::constants::codata::mu_0;

            const double mu_0_val = mu_0.value().value();

            // x(0) -> intitial mu_r
            // x(1) -> saturation polarization

            fvec = 2. * x(1) / Pi()
               * atan(
                  Pi() / 2. / x(1) * mu_0_val * (x(0) - 1)
                  * bi_h_values_.col(0)) - bi_h_values_.col(1);

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
            using Eigen::atan;
            using boost::units::si::constants::codata::mu_0;

            const double mu_0_val = mu_0.value().value();

            const auto help1 = Pi() * mu_0_val * (x(0) - 1.) / x(1);

            fjac.col(0) = 2. * mu_0_val * bi_h_values_.col(0)
               / (Pow<2>(help1 * bi_h_values_.col(0)) + 1.);

            fjac.col(1) = -2. * mu_0_val * (x(0) - 1.) * bi_h_values_.col(0)
               / (x(1) * (Pow<2>(help1 * bi_h_values_.col(0)) + 1.))
               + 2. * atan(help1 * bi_h_values_.col(0)) / Pi();

            return 0;
         }
      };
   }

   //! Use the arctan-function to model the B-H-curve
   /*!
    * This class uses the arctan-function to approximate the B-H-curve
    * as follows:
    *
    * \f{align*}{
    * B(H)&=\mu_{0}H + \frac{2*J_{\mathrm{s}}}{\pi}\arctan\frac{\pi(\mu_{\mathrm{r}}-1)\mu_{0}H}{2J_{\mathrm{s}}}  \\
    * \intertext{with}
    * \mu_{\mathrm{r}} &= \text{initial relative permeability} \\
    * J_{\mathrm{s}} &= \text{saturation polarization}
    * \f}
    *
    * This model has the special property to result in a monotonic
    * \f$\nu\f$-B-curve so that the newton-procedure smoothly converges.
    */
   struct ArctanReluctivity:
      ScalarReluctivity, detail::MakeShared<ArctanReluctivity>
   {
   private:
      template<typename T>
      using AddConstAndRef = detail::AddConstAndRef<T>;

   public:
      //! Initialize a ArctanReluctivity
      /*!
       * @param initial_mu_r is the initial relative permeabiliy
       * @param saturation_polarization is the saturation polarization
       *
       * Use a curve fitting procedure to find suitable values for
       * these parameters.
       */
      ArctanReluctivity(double initial_mu_r,
                        MagneticPolarizationReal saturation_polarization) :
         initial_mu_r_(initial_mu_r), saturation_polarization_(
            saturation_polarization)
      {
      }

      //! Copy construct
      /*!
       * @param source is the object to copy from
       */
      ArctanReluctivity(const ArctanReluctivity& source) :
         initial_mu_r_(source.initial_mu_r_), saturation_polarization_(
            source.saturation_polarization_)
      {
      }

      //! Assignment
      ArctanReluctivity& operator = (const ArctanReluctivity& source) {
         this->initial_mu_r_ = source.initial_mu_r_;
         this->saturation_polarization_ = source.saturation_polarization_;

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
      ArctanReluctivity(HIterator, AddConstAndRef<HIterator>,
                        BiIterator);

      //! Calculate the Reluctivity and its derivation for a given flux density value
      /*!
       * @param b is the magnetic flux density for which the
       * reluctivity and its derivation shall be computed.
       *
       * @returns a pair of the reluctivity and its derivation.
       */
      virtual ReluctivityReluctivityDerivedPair operator()(
         MagneticFluxDensityReal b) const override;

      bool ResultCacheRecommended() const override;

      //! Retrieve the model's initial mu_r
      /*!
       * Especially useful after initializing this class with data points.
       *
       * @return the model's initial mu_r
       */
      double initial_mu_r() const
      {
         return initial_mu_r_;
      }

      //! Retrieve the model's saturation polarization
      /*!
       * Especially useful after initializing this class with data points.
       *
       * @return the model's saturation polarization
       */
      MagneticPolarizationReal saturation_polarization() const
      {
         return saturation_polarization_;
      }

   protected:
      //! Do the initialization
      template <typename HIterator, typename BiIterator>
      void Init(HIterator, AddConstAndRef<HIterator>,
                BiIterator);

   private:
      using DoubleTuple = std::tuple<double, double>;

      //! Helper function to solve the B-H-curve-equation for nu
      /*!
       * Since the arctan-magnetization-curve can not be solved
       * analyticaly, it is solved numericaly using the newton
       * procedure.
       *
       * @param b is the magnetic flux density under consideration
       * @param nu is the reluctivity in the current newton step
       *
       * @return function value and derivative
       */
      DoubleTuple NewtonIteration(MagneticFluxDensityReal, double) const;

   private:
      double initial_mu_r_;
      MagneticPolarizationReal saturation_polarization_;
      mutable std::map<MagneticFluxDensityReal,
                       MagneticReluctivity> value_cache_;
   };

   //--------------------------------------------------------------
   //
   //   Function implementations
   //
   //--------------------------------------------------------------

   template<typename HIterator, typename BiIterator>
   ArctanReluctivity::ArctanReluctivity(
      HIterator h_begin, AddConstAndRef<HIterator> h_end, BiIterator bi_begin)
   {
      Init(h_begin, h_end, bi_begin);
   }

   template <typename HIterator, typename BiIterator>
   inline void ArctanReluctivity::Init(
      HIterator h_begin, AddConstAndRef<HIterator> h_end,
      BiIterator bi_begin)
   {
      using boost::units::si::tesla;

      Eigen::VectorXd x(2);
      // First try
      x << 5e3, 1.8;
      detail::LMFunctorArctan lmfunc(h_begin, h_end, bi_begin);
      Eigen::LevenbergMarquardt<detail::LMFunctorArctan> lm(lmfunc);
      const auto info = lm.minimize(x);
      Eigen::CheckLevenbergMarquardtStatus(info);

      initial_mu_r_ = x(0);
      saturation_polarization_ = x(1) * tesla;
   }
}

#endif /* ARCTAN_RELUCTIVITY_HPP_DvG5eh2X_ */

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
