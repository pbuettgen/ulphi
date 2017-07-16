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

#include "grid_assembly.hpp"

namespace PACKAGE_NS {
   auto GridAssembly::AddGrid(ConstGridPointer grid) -> GridIndex
   {
      grids_.push_back(grid);
      return grids_.size() - 1;
   }

   auto GridAssembly::NumberOfVerticesBrutto() const -> VertexIndex
   {
      VertexIndex total_number_of_vertices = 0;

      for (auto iter = grids_.begin(), end = grids_.end(); iter != end;
           ++iter) {
         total_number_of_vertices += (*iter)->NumberOfVerticesBrutto();
      }

      return total_number_of_vertices;
   }

   auto GridAssembly::NumberOfEdgesNetto() const -> VertexIndex
   {
      VertexIndex total_number_of_edges = 0;

      for (auto iter = grids_.begin(), end = grids_.end(); iter != end;
           ++iter) {
         total_number_of_edges += (*iter)->NumberOfEdgesNetto();
      }

      return total_number_of_edges;
   }

   auto GridAssembly::DoSetFirstVertexIndex(
      VertexIndex vertex_index) const -> VertexIndex
   {
      VertexIndex i = vertex_index;

      for (auto iter = grids_.begin(), end = grids_.end(); iter != end;
           ++iter) {
         i = (*iter)->SetFirstVertexIndex(i);
      }

      return i;
   }
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
