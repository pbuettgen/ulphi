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
 *  \brief A rectangle
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef RECTANGLE_HPP_7AfBpwLc_
#define RECTANGLE_HPP_7AfBpwLc_

#include <tuple>

#include <boost/format.hpp>

#include "config.h"

#include "exception.hpp"
#include "shape.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   //! Rectangle
   /*!
    * Rectangle means a shape having edges which are parallel to the
    * coordinate system's axis.
    *
    * \tparam CoordinateSystem is the coordinate system type
    */
   template<typename CoordinateSystem>
   struct Rectangle :
      Shape<CoordinateSystem>,
      detail::MakeShared<Rectangle<CoordinateSystem> >
   {
   private:
      using Point = typename CoordinateSystem::Point;
      using FirstCoordinate = typename CoordinateSystem::FirstCoordinate;
      using SecondCoordinate = typename CoordinateSystem::SecondCoordinate;

   public:
      //! Initialization
      /*!
       * @param lower_left is the rectangle's lower left corner
       * @param upper_right is the rectangle's upper right corner
       * @throw InvalidArgument in the case of an invalid definition
       */
      Rectangle(const Point& lower_left, const Point& upper_right) :
         lower_left_(lower_left), upper_right_(upper_right)
      {
         CheckValid();
      }

      //! Initialization
      /*!
       * @param lower_left is the rectangle's lower left corner
       *
       * @param xsize is the rectangle's dimension in the coordinate
       * system's first axis' direction.
       *
       * @param ysize is the rectangle's dimension in the coordinate
       * system's second axis' direction.
       */
      Rectangle(const Point& lower_left,
                FirstCoordinate xsize, SecondCoordinate ysize)
         : lower_left_(lower_left), upper_right_(
            std::make_tuple(std::get<0>(lower_left) + xsize,
                            std::get<1>(lower_left) + ysize))
      {
      }

      //! Does a given point lie within this rectangle?
      bool Within(const Point& point) const override
      {
         using std::get;
         return (get<0>(lower_left_) <= get<0>(point))
            && (get<0>(upper_right_) >= get<0>(point))
            && (get<1>(lower_left_) <= get<1>(point))
            && (get<1>(upper_right_) >= get<1>(point));
      }

      Area ComputeArea() const override
      {
         return CoordinateSystem::SurfaceElement(lower_left_, upper_right_);
      }

   protected:
      void CheckValid() const;

   private:
      Point lower_left_, upper_right_;
   };

   //--------------------------------------------------------------
   //
   //   Function implementations
   //
   //--------------------------------------------------------------

   template<typename CoordinateSystem>
   inline void Rectangle<CoordinateSystem>::CheckValid() const
   {
      using std::get;

      if ((get<0>(lower_left_) >= get<0>(upper_right_))
          or (get<1>(lower_left_) >= get<1>(upper_right_))) {
#ifdef __GNUC__
#define FUNCTION_NAME __PRETTY_FUNCTION__
#else
#define FUNCTION_NAME (typeid(*this).name())
#endif
         boost::format error_message("%1%: Invalid rectangle definition!");
         throw InvalidArgument((error_message % FUNCTION_NAME).str());
#undef FUNCTION_NAME
      }
   }

#ifndef INSTANTIATE_
   extern template struct Rectangle<CartesianSystem>;
   extern template struct Rectangle<PolarSystem>;
#endif
}

#endif // RECTANGLE_HPP_7AfBpwLc_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
