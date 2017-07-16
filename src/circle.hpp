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
 *  \brief A circle
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef CIRCLE_HPP_f6qpSpEk_
#define CIRCLE_HPP_f6qpSpEk_

#include "config.h"

#include <boost/math/constants/constants.hpp>
#include <boost/units/cmath.hpp>

#include "coordinate_systems.hpp"
#include "shape.hpp"
#include "pow.hpp"
#include "types.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   //! A circle
   /*!
    * @tparam CoordinateSystem is the coordinate system the circle is
    * defined in.
    */
   template <typename CoordinateSystem>
   struct Circle :
      Shape<CoordinateSystem>,
      detail::MakeShared<Circle<CoordinateSystem> >
   {
   private:
      //! Point type
      using Point = typename CoordinateSystem::Point;

   public:
      //! Initialization
      /*!
       * @param centre is the circle's center point
       * @param radius is the circle's radius
       */
      Circle(const Point& centre, Length radius) :
            centre_(centre), radius_(abs(radius))
      {
      }

      //! Does a given point lie within this circle?
      /*!
       * A point onto the circle's boundary also counts as a point
       * lying within the circle.
       *
       * @param point is the point to check
       * @return whether the passed point lies within the circle.
       */
      bool Within(const Point& point) const override
      {
         const auto distance = CoordinateSystem::EuclideanDistance(
	    centre_, point);
         return distance <= radius_;
      }

      Area ComputeArea() const override
      {
	 return boost::math::double_constants::pi * Pow<2>(radius_);
      }

   private:
      Point centre_;
      Length radius_;
   };

#ifndef INSTANTIATE_
   extern template struct Circle<CartesianSystem>;
   extern template struct Circle<PolarSystem>;
#endif
}

#endif // CIRCLE_HPP_f6qpSpEk_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
