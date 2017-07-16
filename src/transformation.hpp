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
 *  \brief Transform a point
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef TRANSFORMATION_HPP_3D695Ran_
#define TRANSFORMATION_HPP_3D695Ran_

#include <memory>

#include "config.h"

#include "coordinate_systems.hpp"
#include "types.hpp"

namespace PACKAGE_NS {
   //! Base class for all transformations
   struct Transformation {
      using Point = CartesianSystem::Point;

   protected:
      using TransformationPointer = std::unique_ptr<Transformation>;

   public:
      //! Destruction
      virtual ~Transformation();

      //! Apply the transformation
      /*!
       * @param point is the point to transform.
       */
      virtual void Apply(Point&) const =0;

      virtual TransformationPointer Clone() const =0;
   };

   //! Shift
   struct Shift: Transformation {
      using Vector = CartesianSystem::Vector;

      //! Initialization
      /*!
       * @param shift_vector shift vector
       */
      /*
        Shift(const CartesianSystem::Vector& shift_vector) :
        shift_vector_(
        MapVector(
        const_cast<CartesianSystem::DataType*>(&(std::get<0>(
        shift_vector).value())), CartesianSystem::DIMENSION, 1))
        {
        }
      */

      Shift(const Vector& shift_vector) :
         shift_vector_(shift_vector)
      {
      }

      //! Destruction
      virtual ~Shift();

      void Apply(Point&) const override;

      TransformationPointer Clone() const override;

   private:
      Vector shift_vector_;
   };

   //! Rotation
   struct Rotate: Transformation {
      using Float = Angle::value_type;
      //! Initialization
      /*!
       * @param angle rotation angle
       */
      Rotate(Angle angle) :
         cos_alpha_(cos(angle)), sin_alpha_(sin(angle))
         /*angle_(angle)*/
      {
      }

      //! Destruction
      virtual ~Rotate();

      void Apply(Point&) const override;

      TransformationPointer Clone() const override;

   private:
      Float cos_alpha_;
      Float sin_alpha_;
//      Angle angle_;
   };
}

#endif // TRANSFORMATION_HPP_3D695Ran_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
