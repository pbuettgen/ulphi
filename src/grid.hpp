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
 *  \brief Interface for all grid types
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef GRID_HPP_d296FE07_
#define GRID_HPP_d296FE07_

#include "config.h"

#include "loki_visitor_extra.hpp"
#include "types.hpp"

namespace PACKAGE_NS {
   //! Interface class for all grids
   struct Grid : Loki::BaseVisitable<void, Loki::DefaultCatchAll, true>
   {
   protected:
      //! Initialize a grid
      Grid() : first_vertex_index_(0)
      {
      }

   public:
      //! Destruct a grid
      virtual ~Grid()
      {
      }

      //! Compute the number of vertices within this grid
      /*!
       * Also count vertices hidden by a boundary coupling
       *
       * @return total number of vertices
       */
      virtual VertexIndex NumberOfVerticesBrutto() const =0;

      //! Total number of edges
      /*!
       * Edges hidden by grid couplings are not counted.
       *
       * @return total number of edges
       */
      virtual VertexIndex NumberOfEdgesNetto() const =0;

      //! Set the first vertex number this grid uses
      /*!
       * @param vertex_number is the first vertex number used by this
       * grid.
       *
       * @return the next free global vertex index
       */
      VertexIndex SetFirstVertexIndex(VertexIndex vertex_index) const
      {
         first_vertex_index_ = vertex_index;
         return DoSetFirstVertexIndex(vertex_index);
      }

      VertexIndex SetFirstVertexIndex() const {
      	 return SetFirstVertexIndex(0);
      }

   protected:
      //! Assign first vertex index hook
      virtual VertexIndex DoSetFirstVertexIndex(
         VertexIndex vertex_index) const =0;

      //! Convert a local vertex index to a global one
      /*!
       * @param local_vertex_index is the local vertex index
       * @return local vertex index + first vertex index
       */
      VertexIndex ConvertLocal2Global(VertexIndex local_vertex_index) const
      {
         return first_vertex_index_ + local_vertex_index;
      }

      //! Convert a global index to a local one
      VertexIndex ConvertGlobal2Local(VertexIndex global_vertex_index) const
      {
         return global_vertex_index - first_vertex_index_;
      }

   private:
      mutable VertexIndex first_vertex_index_;
   };
}

#endif // GRID_HPP_d296FE07_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */

