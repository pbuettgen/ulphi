// -*- coding: utf-8 -*-
/*
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and library for
 * computing electromagnetic fields.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*!
 *  \file
 *
 *  \brief Shape base class
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef SHAPE_HPP_9LLwxwk0_
#define SHAPE_HPP_9LLwxwk0_

#include "config.h"

namespace PACKAGE_NS {
   //! Shape base class
   /*!
    * \tparam CoordinateSystem is the coordinate system type
    */
   template<typename CoordinateSystem>
   struct Shape {
      using Point = typename CoordinateSystem::Point;

      //! Destruction
      virtual ~Shape();

      //! Does a given point lie within this shape?
      virtual bool Within(const Point&) const =0;

      //! Compute a shape's area
      virtual Area ComputeArea() const =0;
   };

   template <typename CoordinateSystem>
   using ConstShapePointer = std::shared_ptr<const Shape<CoordinateSystem> >;

   //--------------------------------------------------------------
   //
   //   Function implementations
   //
   //--------------------------------------------------------------

   template<typename CoordinateSystem>
   Shape<CoordinateSystem>::~Shape()
   {
   }
}

#endif // SHAPE_HPP_9LLwxwk0_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */

