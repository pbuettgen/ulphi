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

#include "boost_units_extra.hpp"
#include "constant_reluctivity.hpp"

auto PACKAGE_NS::ConstantReluctivity::operator ()(
   MagneticFluxDensityReal b) const -> ReluctivityReluctivityDerivedPair
{
   using boost::units::si::metre_per_henry_per_tesla;

   return ReluctivityReluctivityDerivedPair(nu_, .0*metre_per_henry_per_tesla);
}

bool PACKAGE_NS::ConstantReluctivity::ResultCacheRecommended() const
{
   return false;
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
