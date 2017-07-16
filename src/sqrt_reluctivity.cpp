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

#include <cmath>
#include <functional>
#include <limits>

#include <boost/format.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/units/cmath.hpp>

#include "pow.hpp"
#include "exception.hpp"
#include "sqrt_reluctivity.hpp"

namespace PACKAGE_NS {
   ReluctivityReluctivityDerivedPair SqrtReluctivity::CalculateWithNewton(
      MagneticFluxDensityReal b) const
   {
      using std::placeholders::_1;
      using boost::units::si::metre_per_henry;
      using boost::units::si::tesla;
      using boost::units::si::constants::codata::mu_0;
      using boost::math::tools::newton_raphson_iterate;

      b = abs(b);

      // Avoid division by zero
      if (.0 * tesla == b)
         b = 1e-4 * tesla;

      // Calculate reluctivity corresponding to b
      const auto cache_iter = value_cache_.lower_bound(b);
      const bool cache_success = (value_cache_.end() != cache_iter);
      const double reluctivity_estimate = (
         cache_success ? cache_iter->second : 1. / mu_0 / initial_mu_r_)
         / (1. * metre_per_henry);
      constexpr std::uintmax_t kIterLimit = 64;
      std::uintmax_t max_iterations = kIterLimit;
      const MagneticReluctivity nu = metre_per_henry
         * newton_raphson_iterate(
            std::bind(&SqrtReluctivity::NewtonIteration, this, b, _1),
            reluctivity_estimate, .0, 20.*reluctivity_estimate,
            std::numeric_limits<float>::digits10, max_iterations);

      if (max_iterations >= kIterLimit)
      {
         boost::format error_message(
            "SqrtReluctivity: No convergence in reluctivity calculation after %1% iterations!");
         throw RuntimeError((error_message % max_iterations).str());
      }

      if ((not cache_success)
          or (cache_success and (
                 (b / cache_iter->first < .9)
                 or (nu/cache_iter->second < .9))))
         value_cache_.emplace_hint(cache_iter, b, nu);

      const auto nu_diff =
         (-1. * mu_0 * nu
          - saturation_polarization_
          * (mu_0 * nu * (initial_mu_r_ - 1.)
             / saturation_polarization_
             - (-2. * mu_0 * nu * (-a_ + 1.)
                * (initial_mu_r_ - 1.)
                / saturation_polarization_
                + mu_0 * nu * (initial_mu_r_ - 1.)
                * (mu_0 * b * nu
                   * (initial_mu_r_ - 1.)
                   / saturation_polarization_ + 1.)
                / saturation_polarization_)
             / sqrt(
                -4. * mu_0 * b * nu * (-a_ + 1.)
                * (initial_mu_r_ - 1.)
                / saturation_polarization_
                + Pow<2>(
                   mu_0 * b * nu
                   * (initial_mu_r_
                      - 1.)
                   / saturation_polarization_
                   + 1.)))
          / (-2. * a_ + 2.) + 1.)
         / (mu_0 * b
            + saturation_polarization_
            * (mu_0 * b * (initial_mu_r_ - 1.)
               / saturation_polarization_
               - (-2. * mu_0 * b * (-a_ + 1.)
                  * (initial_mu_r_ - 1.)
                  / saturation_polarization_
                  + mu_0 * b * (initial_mu_r_ - 1.)
                  * (mu_0 * b * nu
                     * (initial_mu_r_ - 1.)
                     / saturation_polarization_
                     + 1.)
                  / saturation_polarization_)
               / sqrt(
                  -4. * mu_0 * b * nu * (-a_ + 1.)
                  * (initial_mu_r_ - 1.)
                  / saturation_polarization_
                  + Pow<2>(
                     mu_0 * b * nu
                     * (initial_mu_r_
                        - 1.)
                     / saturation_polarization_
                     + 1.)))
            / (-2. * a_ + 2.));

      return std::make_tuple(nu, nu_diff);
   }

   ReluctivityReluctivityDerivedPair SqrtReluctivity::CalculateDirectly(
      MagneticFluxDensityReal b) const
   {
      using boost::units::si::metre_per_henry;
      using boost::units::si::tesla;
      using boost::units::si::constants::codata::mu_0;

      // Avoid division by zero
      if (.0 * tesla == b)
         b = 1e-4 * tesla;

      const MagneticReluctivity nu = .5
         * (mu_0
            * (2. * a_ * b - b * initial_mu_r_ - b
               + initial_mu_r_ * saturation_polarization_)
            - sqrt(
               Pow<2>(mu_0)
               * (4. * a_ * b * initial_mu_r_
                  * saturation_polarization_
                  - 4. * a_ * b * saturation_polarization_
                  + Pow<2>(b * initial_mu_r_)
                  - 2. * Pow<2>(b) * initial_mu_r_ + Pow<2>(b)
                  - 2. * b * Pow<2>(initial_mu_r_)
                  * saturation_polarization_
                  + 2. * b * initial_mu_r_
                  * saturation_polarization_
                  + Pow<2>(
                     initial_mu_r_
                     * saturation_polarization_))))
         / (b * Pow<2>(mu_0) * (a_ - initial_mu_r_));

      const auto nu_diff =
         .5
         * (mu_0 * (2. * a_ - initial_mu_r_ - 1)
            - .5
            * sqrt(
               Pow<2>(mu_0)
               * (4. * a_ * b * initial_mu_r_
                  * saturation_polarization_
                  - 4. * a_ * b
                  * saturation_polarization_
                  + Pow<2>(b * initial_mu_r_)
                  - 2. * Pow<2>(b) * initial_mu_r_
                  + Pow<2>(b)
                  - 2. * b * Pow<2>(initial_mu_r_)
                  * saturation_polarization_
                  + 2. * b * initial_mu_r_
                  * saturation_polarization_
                  + Pow<2>(
                     initial_mu_r_
                     * saturation_polarization_)))
            * (4. * a_ * initial_mu_r_
               * saturation_polarization_
               - 4. * a_ * saturation_polarization_
               + 2. * b * Pow<2>(initial_mu_r_)
               - 4. * b * initial_mu_r_ + 2. * b
               - 2. * Pow<2>(initial_mu_r_)
               * saturation_polarization_
               + 2. * initial_mu_r_
               * saturation_polarization_)
            / (4. * a_ * b * initial_mu_r_
               * saturation_polarization_
               - 4. * a_ * b * saturation_polarization_
               + Pow<2>(b * initial_mu_r_)
               - 2. * Pow<2>(b) * initial_mu_r_ + Pow<2>(b)
               - 2. * b * Pow<2>(initial_mu_r_)
               * saturation_polarization_
               + 2. * b * initial_mu_r_
               * saturation_polarization_
               + Pow<2>(
                  initial_mu_r_
                  * saturation_polarization_)))
         / (b * Pow<2>(mu_0) * (a_ - initial_mu_r_))
         - .5
         * (mu_0
            * (2. * a_ * b - b * initial_mu_r_ - b
               + initial_mu_r_ * saturation_polarization_)
            - sqrt(
               Pow<2>(mu_0)
               * (4. * a_ * b * initial_mu_r_
                  * saturation_polarization_
                  - 4. * a_ * b
                  * saturation_polarization_
                  + Pow<2>(b * initial_mu_r_)
                  - 2. * Pow<2>(b) * initial_mu_r_
                  + Pow<2>(b)
                  - 2. * b * Pow<2>(initial_mu_r_)
                  * saturation_polarization_
                  + 2. * b * initial_mu_r_
                  * saturation_polarization_
                  + Pow<2>(
                     initial_mu_r_
                     * saturation_polarization_))))
         / (Pow<2>(b * mu_0) * (a_ - initial_mu_r_));

      return std::make_tuple(nu, nu_diff);
   }

   ReluctivityReluctivityDerivedPair SqrtReluctivity::operator ()(
      MagneticFluxDensityReal b) const
   {
      return this->CalculateWithNewton(b);
   }

   bool SqrtReluctivity::ResultCacheRecommended() const
   {
      return true;
   }

   auto SqrtReluctivity::NewtonIteration(MagneticFluxDensityReal b,
                                         double nu_without_unit) const -> DoubleTuple
   {
      using std::sqrt;
      using boost::units::si::metre_per_henry;
      using boost::units::si::tesla;
      using boost::units::si::constants::codata::mu_0;

      const auto nu = nu_without_unit * metre_per_henry;

      const auto help1 = 1. - a_;
      const auto help2 = initial_mu_r_ - 1.;
      const auto help3 = mu_0 * b;
      const auto help4 = help2 * help3;
      const auto h = nu * b;
      const double h_a = mu_0 * h * help2 / saturation_polarization_;

      const auto delta = mu_0 * h - b
         + saturation_polarization_
         * ((h_a + 1 - sqrt(Pow<2>(h_a + 1) - 4 * h_a * help1))
            / (2 * help1));

      const auto diff =
         help3
         + saturation_polarization_
         * (help4 / saturation_polarization_
            - (-2. * help4 * help1 / saturation_polarization_
               + help4
               * (help4 * nu
                  / saturation_polarization_ + 1)
               / saturation_polarization_)
            / sqrt(
               -4. * help1
               * (help4 * nu
                  / saturation_polarization_)
               + Pow<2>(
                  help4 * nu
                  / saturation_polarization_
                  + 1)))
         / (2 * help1);

      return std::make_tuple(delta / tesla, diff / tesla * metre_per_henry);
   }

}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
