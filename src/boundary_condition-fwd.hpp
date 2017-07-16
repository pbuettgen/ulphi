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
 *  \brief Boundary conditions -- forward declarations
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef BOUNDARY_CONDITION_FWD_HPP_q8uYYmkT_
#define BOUNDARY_CONDITION_FWD_HPP_q8uYYmkT_

#include <memory>

#include "config.h"

namespace PACKAGE_NS {
   //! Base class
   template <typename > struct BoundaryCondition;

   //! Dirichlet condition
   template <typename > struct DirichletCondition;

   //! Inhomogenious neumann condition
   template <typename > struct NeumannCondition;

   //! The void boundary condition
   template <typename > struct VoidCondition;

   //! Pointer to a boundary condition object
   template <typename Analysis>
   using ConstBoundaryConditionPointer =
      std::shared_ptr<const BoundaryCondition<Analysis> >;
}

#endif // BOUNDARY_CONDITION_FWD_HPP_q8uYYmkT_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
