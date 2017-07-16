// -*- coding: utf-8 -*-
/*
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and
 * library for computing electromagnetic fields.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "config.h"

#include "boost_python_extra.hpp"

#include "exp_vertex_distribution.hpp"
#include "linear_vertex_distribution.hpp"
#include "types.hpp"

namespace PACKAGE_NS {
   namespace python {
      // -------------------------
      // *** VertexDistributions ***
      // -------------------------

      struct VertexDistributionWrap
         : VertexDistribution, bp::wrapper<VertexDistribution>
      {
         using Base = VertexDistribution;
         using VertexIndex = Base::VertexIndex;

         VertexDistributionWrap(VertexIndex number_of_vertices)
            : Base(number_of_vertices) {}

         double PositionNthVertex(VertexIndex vertex_index) const override
         {
            return this->get_override("PositionNthVertex")();
         }
      };
   }
}

template <typename VertexDistributionType>
using BPClassVertexDist = boost::python::class_<
   VertexDistributionType,
   std::shared_ptr<VertexDistributionType>,
   boost::python::bases<PACKAGE_NS::VertexDistribution>,
   boost::noncopyable>;

BOOST_PYTHON_MODULE(MODULE_NAME)
{
   using namespace PACKAGE_NS;
   using namespace PACKAGE_NS::python;
   namespace bp = boost::python;

   BPClassSharedPtr<VertexDistributionWrap>(
      "VertexDistribution", bp::no_init)
      .def("PositionNthVertex",
           bp::pure_virtual(&VertexDistribution::PositionNthVertex))
      .add_property("number_of_vertices",
                    &VertexDistributionWrap::number_of_vertices,
                    &VertexDistributionWrap::set_number_of_vertices);

   BPClassVertexDist<LinearVertexDistribution>(
      "LinearVertexDistribution", bp::init<VertexIndex>())
      .def("PositionNthVertex",
           &LinearVertexDistribution::PositionNthVertex);

   BPClassVertexDist<ExpVertexDistribution>(
      "ExpVertexDistribution", bp::init<VertexIndex>())
      .def("PositionNthVertex",
           &ExpVertexDistribution::PositionNthVertex);

   bp::implicitly_convertible<
      std::shared_ptr<LinearVertexDistribution>,
      ConstVertexDistributionPointer>();
   bp::implicitly_convertible<
      std::shared_ptr<ExpVertexDistribution>,
      ConstVertexDistributionPointer>();
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
