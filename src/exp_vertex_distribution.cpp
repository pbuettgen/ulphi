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

#include <cmath>
#include "exp_vertex_distribution.hpp"

namespace PACKAGE_NS {
   double ExpVertexDistribution::PositionNthVertex(
      VertexIndex vertex_index) const
   {
      CheckValid(vertex_index);

      return
         reverse_ ?
         (this->number_of_vertices() - 1 == vertex_index ?
          1. : 1. - std::pow(base_, -static_cast<double>(vertex_index))) :
         (0 == vertex_index ?
          0. :
          std::pow(base_,
                   static_cast<double>(vertex_index) + 1.
                   - this->number_of_vertices()));
   }
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
