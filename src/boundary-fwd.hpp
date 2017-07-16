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
 *  \brief A GridPatch's boundary -- forward declarations
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef BOUNDARY_FWD_HPP_JbT4DfOH_
#define BOUNDARY_FWD_HPP_JbT4DfOH_

#include "config.h"

namespace PACKAGE_NS {
   //! Boundary base class
   template<typename > struct Boundary;

   //! East boundary
   template<typename > struct EastBoundary;

   //! North boundary
   template<typename > struct NorthBoundary;

   //! West boundary
   template<typename > struct WestBoundary;

   //! South boundary
   template<typename > struct SouthBoundary;
}

#endif // BOUNDARY_FWD_HPP_JbT4DfOH_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
