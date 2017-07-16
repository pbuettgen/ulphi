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

#include "boundary.hpp"

namespace PACKAGE_NS {
   template struct EastBoundary<AnalysisType>;
   template struct NorthBoundary<AnalysisType>;
   template struct WestBoundary<AnalysisType>;
   template struct SouthBoundary<AnalysisType>;
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
