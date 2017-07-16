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

#include <functional>
#include <limits>
#include <cstdint>

#include <boost/format.hpp>
#include <boost/units/cmath.hpp>

#include "arctan_reluctivity.hpp"
#include "exception.hpp"
#include "pow.hpp"
#include "types.hpp"

namespace PACKAGE_NS {
   ReluctivityReluctivityDerivedPair ArctanReluctivity::operator ()(
      MagneticFluxDensityReal b) const
   {
      using std::placeholders::_1;
      using boost::units::si::metre_per_henry;
      using boost::units::si::tesla;
      using boost::units::si::constants::codata::mu_0;
      using boost::math::tools::newton_raphson_iterate;
      static double kPi = boost::math::double_constants::pi;

      b = boost::units::abs(b);

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
            std::bind(&ArctanReluctivity::NewtonIteration, this, b, _1),
            reluctivity_estimate, 0., 20. * reluctivity_estimate,
            std::numeric_limits<float>::digits10, max_iterations);

      if (max_iterations >= kIterLimit)
      {
         boost::format error_message(
            "ArctanReluctivity: No convergence in reluctivity calculation after %1% iterations!");
         throw RuntimeError((error_message % max_iterations).str());
      }

      if ((not cache_success)
          or (cache_success and (
                 (b / cache_iter->first < .9)
                 or (nu/cache_iter->second < .9))))
         value_cache_.emplace_hint(cache_iter, b, nu);

      // Calculate the derivation with respect to b
      const auto helper1 = initial_mu_r_ - 1.;
      const auto helper2 = mu_0
         * (1
            + helper1
            / (1.
               + Pow<2>(
                  kPi / 2 * helper1 * (mu_0 * nu)
                  * (b / saturation_polarization_))));
      const MagneticReluctivityDerived nu_derived = (1. / helper2 - nu) / b;

      return std::make_pair(nu, nu_derived);
   }

   bool ArctanReluctivity::ResultCacheRecommended() const
   {
      return true;
   }

   auto ArctanReluctivity::NewtonIteration(
      MagneticFluxDensityReal b, double nu_without_unit) const -> DoubleTuple
   {
      using boost::units::si::metre_per_henry;
      using boost::units::si::tesla;
      using boost::units::si::constants::codata::mu_0;
      static double kPi = boost::math::double_constants::pi;

      const MagneticReluctivity nu = nu_without_unit * metre_per_henry;
      const double mur_1 = initial_mu_r_ - 1.;
      const auto atan_arg_fac
         = kPi/2./saturation_polarization_*mur_1*mu_0;

      const MagneticFluxDensityReal delta =
         b * ((mu_0*nu) - 1.) + 2./kPi * saturation_polarization_
         * atan(static_cast<double>(atan_arg_fac*nu*b));

      const auto diff = mu_0 * b
         * (1. + mur_1/(1.+Pow<2>(atan_arg_fac*nu*b)));

      return std::make_tuple(delta / tesla, diff / tesla * metre_per_henry);
   }
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
