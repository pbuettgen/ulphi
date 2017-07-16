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
 *  \brief Coordinate system definitions
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef COORDINATE_SYSTEMS_HPP_cRIhFmA6_
#define COORDINATE_SYSTEMS_HPP_cRIhFmA6_

#include <utility>
#include <tuple>

#include "config.h"

#include <boost/units/systems/si.hpp>
#include <boost/units/cmath.hpp>

#include "exception.hpp"
#include "pow.hpp"
#include "types.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {

   //! Is an edge parallel or antiparallel?
   enum struct Parallelity: short { kAntiparallel=-1, kParallel=1 };

   //! Axes
   enum struct Axis: unsigned short { kFirst=0, kSecond=1 };

   //! Flip an axis
   inline Axis Flip(Axis axis) {
      return Axis::kFirst == axis ? Axis::kSecond : Axis::kFirst;
   }

   //! Compute midpoint
   /*!
    *  Walk only parallel to the coordinate system's axes.
    *
    *  @param p1 first point
    *
    *  @param p2 second point
    *
    *  @return computed midpoint
    */
   template <typename Point>
   Point ManhattenMidpoint(Point p1, Point p2)
   {
      using std::get;

      return std::make_tuple(
         (get<0>(p1)+get<0>(p2))/2., (get<1>(p1)+get<1>(p2))/2.);
   }

   //! Define a cartesian coordinate system's properties
   /*!
    *  The properties for doing computations in cartesian coordinates
    *  are provided.
    */
   struct CartesianSystem {
      //! Dimension of this coordinate system
      static constexpr std::size_t DIMENSION = 2;

      //! data type for (almost) all computations
      using DataType = double;

      //! Type used for processing the x coordinate.
      using FirstCoordinate = boost::units::quantity<
         boost::units::si::length, DataType>;

      //! Type used for processing the y coordinate.
      using SecondCoordinate = boost::units::quantity<
         boost::units::si::length, DataType>;

      //! Type describing a point in this system
      using Point = std::tuple<FirstCoordinate, SecondCoordinate>;

      //! Position vector type
      using PositionVector = Point;

      //! Vector type
      using Vector = Point;

      using VectorNormalized = std::tuple<DataType, DataType>;

      //! Convert to cartesian coordinates
      /*!
       * @param p is the point to convert
       * @return p converted to cartesian coordinates
       */
      static Point ToCartesic(const Point& p)
      {
         return p;
      }

      //! Compute a surface element's area
      /*!
       *  @param lower_left_corner lower left corner
       *
       *  @param upper_right_corner upper right corner
       *
       *  @returns area
       *
       *  @throw InvalidArgument if an edge length is less than zero
       *  (not NDEBUG only)
       */
      static Area SurfaceElement(const Point& lower_left_corner,
                                 const Point& upper_right_corner)
      {
         using std::get;

#ifndef NDEBUG
         if ((get<0>(lower_left_corner) > get<0>(upper_right_corner))
             || (get<1>(lower_left_corner) > get<1>(upper_right_corner)))
            throw RuntimeError(
               PACKAGE_NAME ": Invalid coordinates in cell size computation!");
#endif
         return (get<0>(upper_right_corner)-get<0>(lower_left_corner))
            * (get<1>(upper_right_corner)-get<1>(lower_left_corner));
      }

      //! Compute the euclidean distance between two points
      /*!
       * Euclidean distance means the length of a straight line
       * between the given points.
       *
       * @param p1 is point 1
       * @param p2 is point 2
       * @return the euclidean distance between the points
       */
      static Length EuclideanDistance(const Point& p1, const Point& p2)
      {
         using std::get;

         return sqrt(
            Pow<2>(get<0>(p2) - get<0>(p1))
            + Pow<2>(get<1>(p2) - get<1>(p1)));
      }

      //! Compute the Manhattan distance between two points
      /*!
       * Manhattan distance means to take the Length when only walking
       * parallel to the coordinate axes.  Walk first only parallel to
       * the first axis, then parallel to the second axis.
       *
       * @param p1 is the first point
       * @param p2 is the second point
       * @return the maximum Manhattan distance between the given points
       */
      static Length ManhattenDistanceMax(const Point& p1, const Point& p2)
      {
         using std::get;

         return abs(get<0>(p2) - get<0>(p1)) + abs(get<1>(p2) - get<1>(p1));
      }
   };

   //! Define a polar coordinate system's properties
   /*!
    *  The properties for doing computations in polar coordinates
    *  are provided.
    */
   struct PolarSystem {
      //! Dimension of this coordinate system
      static constexpr std::size_t DIMENSION=2;

      //! Type for (almost) all computations
      using DataType = double;

      //! Type used for processing the r coordinate.
      using FirstCoordinate = boost::units::quantity<
         boost::units::si::length, DataType>;

      //! Type used for processing the phi coordinate.
      using SecondCoordinate = boost::units::quantity<
         boost::units::si::plane_angle, DataType>;

      //! Vector type
      using Vector = std::tuple<FirstCoordinate, SecondCoordinate>;

      //! Type describing a point in this system
      using PositionVector = Vector;

      //! Position vector type
      using Point = PositionVector;

      //! Convert to cartesian coordinates
      /*!
       * @param p is the point to convert
       * @return p converted to cartesian coordinates
       */
      static CartesianSystem::Point ToCartesic(const Point& p)
      {
         using std::get;
         return std::make_tuple(get<0>(p) * cos(get<1>(p)),
                                get<0>(p) * sin(get<1>(p)));
      }

      //! Compute a surface element's area
      /*!
       *  @param lower_left_corner lower left corner
       *
       *  @param upper_right_corner upper right corner
       *
       *  @returns area
       *
       *  @throw InvalidArgument if an edge length is less than zero
       *  (not NDEBUG only)
       */
      static Area SurfaceElement(const Point& lower_left_corner,
                                 const Point& upper_right_corner)
      {
         using std::get;
         using boost::units::si::radians;

#ifndef NDEBUG
         if ((get<0>(lower_left_corner) > get<0>(upper_right_corner))
             || (get<1>(lower_left_corner) > get<1>(upper_right_corner)))
            throw InvalidArgument(
               PACKAGE_NAME ": Invalid coordinates in cell size computation!");
#endif
         return .5 * (Pow<2>(get<0>(upper_right_corner))
                      -Pow<2>(get<0>(lower_left_corner)))
            * ((get<1>(upper_right_corner)
                -get<1>(lower_left_corner))/radians);
      }

      //! Compute the euclidean distance between two points
      /*!
       * Euclidean distance means the length of a straight line
       * between the given points.
       *
       * @param p1 is point 1
       * @param p2 is point 2
       * @return the euclidean distance between the points
       */
      static Length EuclideanDistance(const Point& p1, const Point& p2)
      {
         using std::get;

         return sqrt(
            Pow<2>(get<0>(p1)) + Pow<2>(get<0>(p2))
            - (get<0>(p1) * get<0>(p2))
            * (2. * cos(get<1>(p2) - get<1>(p1))));
      }

      //! Compute the Manhattan distance between two points
      /*!
       * Manhattan distance means to take the Length when only walking
       * parallel to the coordinate axes.  Walk first only parallel to
       * the first axis, then parallel to the second axis.
       *
       * @param p1 is the first point
       * @param p2 is the second point
       * @return the maximum Manhattan distance between the given points
       */
      static Length ManhattenDistanceMax(const Point& p1, const Point&p2)
      {
         using std::get;
         using boost::units::si::radians;

         const Length l1 = abs(get<1>(p2) - get<1>(p1)) / radians
            * std::max(get<0>(p1), get<0>(p2));
         const Length l2 = abs(get<0>(p2) - get<0>(p1));

         return l1 + l2;
      }
   };

   inline Area ScalarProduct(
      const CartesianSystem::Vector& p1, const CartesianSystem::Vector& p2)
   {
      using std::get;

      return get<0>(p1)*get<0>(p2)+get<1>(p1)*get<1>(p2);
   }

   inline Area VectorProduct(
      const CartesianSystem::Vector& p1, const CartesianSystem::Vector& p2)
   {
      using std::get;

      return get<0>(p1)*get<1>(p2)-get<1>(p1)*get<0>(p2);
   }

   inline CartesianSystem::DataType VectorProduct(
      const CartesianSystem::VectorNormalized& v1,
      const CartesianSystem::VectorNormalized& v2)
   {
      using std::get;

      return get<0>(v1)*get<1>(v2)-get<1>(v1)*get<0>(v2);
   }

   inline bool IsStraightLine(const CartesianSystem::Point& p1,
                              const CartesianSystem::Point& p2,
                              const CartesianSystem::Point& p3)
   {
      using std::get;

      return FloatEqual(in_base_units(
                           get<0>(p1)*(get<1>(p2)-get<1>(p3))
                           +get<0>(p2)*(get<1>(p3)-get<1>(p1))
                           +get<0>(p3)*(get<1>(p1)-get<1>(p2))), .0);
   }

   inline Length ComputeLength(
      const CartesianSystem::Vector& v)
   {
      return sqrt(ScalarProduct(v, v));
   }

   inline CartesianSystem::Vector Substract(
      const CartesianSystem::Vector& v1, const CartesianSystem::Vector& v2)
   {
      using std::get;

      return std::make_tuple(get<0>(v1)-get<0>(v2),
                             get<1>(v1)-get<1>(v2));
   }

   inline CartesianSystem::VectorNormalized Normalize(
      const CartesianSystem::Vector& v1)
   {
      using std::get;

      auto length = ComputeLength(v1);

      return std::make_tuple(get<0>(v1)/length, get<1>(v1)/length);
   }
}

#endif // COORDINATE_SYSTEMS_HPP_cRIhFmA6_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
