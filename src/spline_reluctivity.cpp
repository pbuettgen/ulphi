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

#include <boost/units/cmath.hpp>

#include "spline_reluctivity.hpp"

namespace PACKAGE_NS {
   SplineReluctivity::~SplineReluctivity()
   {
   }

   auto SplineReluctivity::operator ()(
      MagneticFluxDensityReal b) const -> ReluctivityReluctivityDerivedPair
   {
      using boost::units::si::tesla;
      using boost::units::si::metre_per_henry;
      using boost::units::si::metre_per_henry_per_tesla;

      b = abs(b);

      if (b < b_min_) {
	 // linear extrapolation
	 const auto derivatives = spline_->derivatives<1>(0.);
	 return std::make_tuple(
            (derivatives(1) * (b - b_min_) / tesla + derivatives(0))
	    * metre_per_henry, derivatives(1) * metre_per_henry_per_tesla);
      }

      if (b > b_max_) {
	 // linear extrapolation
	 const auto derivatives = spline_->derivatives<1>(1.);
	 return std::make_tuple(
            (derivatives(1) * (b - b_max_) / tesla + derivatives(0))
	    * metre_per_henry, derivatives(1) * metre_per_henry_per_tesla);
      }

      const auto derivatives = spline_->derivatives<1>(UValue(b));

      return std::make_tuple(derivatives(0) * metre_per_henry,
			     derivatives(1) * metre_per_henry_per_tesla);
   }

   bool SplineReluctivity::ResultCacheRecommended() const {
      return true;
   }
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
