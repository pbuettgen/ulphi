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
 *  \brief load definition -- forward declarations
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef LOAD_DEFINITIONFWD_HPP_0T9OlJQ3_
#define LOAD_DEFINITIONFWD_HPP_0T9OlJQ3_

#include <memory>

#include "config.h"

namespace PACKAGE_NS {
   //! Load definition
   template <typename > struct LoadDefinition;

   //! Shared pointer to an immutable load definition object
   template <typename Analysis>
   using ConstLoadDefinitionPointer =
      std::shared_ptr<const LoadDefinition<Analysis> >;
}

#endif // LOAD_DEFINITIONFWD_HPP_0T9OlJQ3_
