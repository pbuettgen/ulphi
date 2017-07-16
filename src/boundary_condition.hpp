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
 *  \brief Boundary condition base class
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef BOUNDARY_CONDITION_HPP_Vj0ox7Dm_
#define BOUNDARY_CONDITION_HPP_Vj0ox7Dm_

#include "config.h"

#include "grid_graph-fwd.hpp"
#include "loki_visitor_extra.hpp"
#include "types.hpp"

namespace PACKAGE_NS {

   //! Base class for a BoundaryCondition
   /*!
    * \tparam Analysis is the analysis type
    */
   template<typename Analysis>
   struct BoundaryCondition
      : Loki::BaseVisitable<void, Loki::DefaultCatchAll, true>
   {
      using ReturnType = typename Loki::BaseVisitable<>::ReturnType;

   protected:
      //! Initialization
      explicit BoundaryCondition(short priority)
	 : priority_(priority)
      {
      }

   public:
      //! Destruct a BoundaryCondition
      virtual ~BoundaryCondition()
      {
      }

      //! The priority
      /*!
       * If two boundary conditions are defined (this can happen at
       * corners) the condition with the highest priority wins.
       */
      short priority() const
      {
	 return priority_;
      }

   private:
      short priority_;
   };
}

#endif // BOUNDARY_CONDITION_HPP_Vj0ox7Dm_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
