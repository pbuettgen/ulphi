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
 *  \brief Group several grids together
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef GRID_ASSEMBLY_HPP_9hlYcQZs_
#define GRID_ASSEMBLY_HPP_9hlYcQZs_

#include <initializer_list>
#include <memory>
#include <vector>

#include "config.h"

#include "analysis.hpp"
#include "grid.hpp"
#include "grids-fwd.hpp"
#include "edge_properties.hpp"
#include "vertex.hpp"

namespace PACKAGE_NS {
   //! Grid assembly
   /*!
    * Group several grids together
    *
    *  \tparam Analysis is the analysis type
    */
   struct GridAssembly:
      Grid, detail::MakeShared<GridAssembly>
   {
   private:
      using Base = Grid;

   public:
      using ReturnType = typename Base::ReturnType;

      LOKI_DEFINE_CONST_VISITABLE()

      private:
      template <typename T>
      using AddConstAndRef = typename detail::AddConstAndRef<T>;
      using GridVector = std::vector<ConstGridPointer>;
      using GridIterator = typename GridVector::iterator;
      using GridIndex = typename GridVector::size_type;

   public:
      //! Construct a grid assembly
      GridAssembly()
      {
      }

      //! Initialization
      GridAssembly(std::initializer_list<ConstGridPointer> grids) :
         grids_(grids)
      {
      }

      //! Initialize from a sequence
      /*!
       * @param begin marks the sequence'es start
       *
       * @param end marks the sequence'es end
       */
      template <typename Iter>
      GridAssembly(Iter begin, AddConstAndRef<Iter> end)
         : grids_(begin, end)
      {
      }

      //! Add a grid to this grid assembly
      /*!
       *  \param grid is the grid to add to this assembly
       *
       *  \returns the grid's index within this grid assembly
       */
      GridIndex AddGrid(ConstGridPointer grid);

      //! Get access to all grids
      const GridVector& grids() const {
         return grids_;
      }

      //! Total number of vertices
      /*!
       * Vertices hidden by grid couplings are also counted.
       *
       * @return total number of vertices
       */
      VertexIndex NumberOfVerticesBrutto() const override;

      //! Total number of edges
      /*!
       * Edges hidden by grid couplings are not counted.
       *
       * @return total number of edges
       */
      VertexIndex NumberOfEdgesNetto() const override;

   protected:
      //! Do required extra work when the grid's first vertex index is set
      VertexIndex DoSetFirstVertexIndex(VertexIndex) const override;

   private:
      GridVector grids_;
   };
}

#endif // GRID_ASSEMBLY_HPP_9hlYcQZs_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */

