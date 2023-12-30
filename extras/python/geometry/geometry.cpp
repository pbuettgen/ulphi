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

#include <iostream>

#include "config.h"

#include "boost_python_extra.hpp"

#include "circle.hpp"
#include "rectangle.hpp"

namespace PACKAGE_NS {
   namespace python {
      // -------------
      // *** Shape ***
      // -------------
      template <typename CoordinateSystem>
      struct ShapeWrap
         : Shape<CoordinateSystem>,
         bp::wrapper<Shape<CoordinateSystem> >
      {
         using Point = typename CoordinateSystem::Point;

         bool Within(const Point& p) const override
         {
            bp::override within = this->get_override("Within");

            return within(p);
         }

         Area ComputeArea() const override
         {
            bp::override compute_area = this->get_override("ComputeArea");

            return compute_area();
         }
      };


      template <typename CoordinateSystem>
      void AddClasses()
      {
         namespace bp = boost::python;

         using CSPoint = typename CoordinateSystem::Point;
         using ShapeClass = ShapeWrap<CoordinateSystem>;

         // *** Shape ***
         BPClassSharedPtr<ShapeClass >("Shape", bp::no_init)
            .def("Within", bp::pure_virtual(&ShapeClass::Within))
            .def("ComputeArea", bp::pure_virtual(&ShapeClass::ComputeArea));

         // *** Circle ***
         BPClassSharedPtr<Circle<CoordinateSystem> >("Circle", bp::no_init)
            .def(bp::init<CSPoint, Length>());

         // *** Rectangle ***
         BPClassSharedPtr<Rectangle<CoordinateSystem> >(
            "Rectangle", bp::no_init)
            .def(bp::init<CSPoint, CSPoint>())
            .def(bp::init<CSPoint,
                 typename CoordinateSystem::FirstCoordinate,
                 typename CoordinateSystem::SecondCoordinate>());
      }
   }
}

BOOST_PYTHON_MODULE(MODULE_NAME) {
   using namespace PACKAGE_NS;
   using namespace PACKAGE_NS::python;
   namespace bp = boost::python;
   namespace mpl = boost::mpl;

   AddClasses<CSType>();

   bp::implicitly_convertible<std::shared_ptr<Circle<CSType> >,
                              ConstShapePointer<CSType> >();
   bp::implicitly_convertible<std::shared_ptr<Rectangle<CSType> >,
                              ConstShapePointer<CSType> >();
}
