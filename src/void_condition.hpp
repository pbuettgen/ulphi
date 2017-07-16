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
 *  \brief Boundary with no condition
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef VOID_CONDITION_HPP_Bpnvy5PC_
#define VOID_CONDITION_HPP_Bpnvy5PC_

#include "config.h"

#include "boundary_condition.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   //! Void boundary condition
   /*!
    * This boundary condition is used for all boundaries which have no
    * condition explicitly set.  Effectly it is nothing else then the
    * homogenious neumann boundary condition.
    *
    * @tparam Analysis is the analysis type.
    */
   template <typename Analysis>
   struct VoidCondition :
      BoundaryCondition<Analysis>,
      detail::MakeShared<VoidCondition<Analysis> >
   {
   private:
      using Base = BoundaryCondition<Analysis>;

   public:
      using ReturnType = typename Base::ReturnType;

      LOKI_DEFINE_CONST_VISITABLE()

   private:
      static constexpr short kBasePriority = 0;

   public:
      //! Initialization
      VoidCondition() : Base(kBasePriority)
      {
      }
   };
}

#endif // VOID_CONDITION_HPP_Bpnvy5PC_
