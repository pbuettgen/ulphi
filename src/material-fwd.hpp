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
 *  \brief Describe material properties -- forward declarations
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef MATERIALFWD_HPP_v1YCe06b_
#define MATERIALFWD_HPP_v1YCe06b_

#include <memory>

#include "config.h"

namespace PACKAGE_NS {
   //! Material properties
   struct Material;
   //! Shared pointer to an immutable material object
   using ConstMaterialPointer = std::shared_ptr<const Material>;
}

#endif // MATERIAL_HPP_v1YCe06b_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
