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
 *  \brief Forward declarations for all grids
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef GRIDS_FWD_HPP_3e2Tsj1o_
#define GRIDS_FWD_HPP_3e2Tsj1o_

#include <memory>

#include "config.h"

namespace PACKAGE_NS {
   struct Grid;
   struct GridAssembly;
   template<typename > struct GridPatchBase;
   template<typename, typename > struct GridPatch;

   //! Pointer to a constant grid object
   /*!
    * \tparam Analysis is the analysis type
    *
    * \relates Grid
    */
   using ConstGridPointer = std::shared_ptr<const Grid>;
}

#endif // GRIDS_FWD_HPP_3e2Tsj1o_
