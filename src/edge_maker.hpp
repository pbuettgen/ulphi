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
 *  \brief Call grid edge property objects into life.
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */


#ifndef EDGE_MAKER_HPP_Dengidat_
#define EDGE_MAKER_HPP_Dengidat_

#include "config.h"

#include "grid_graph-fwd.hpp"

namespace PACKAGE_NS {
   namespace fi {
      //! Create a grid edge
      /*!
       *  For polar grids depending on the direction straight or
       *  circular edges are created.
       *
       *  \tparam Analysis is the analysis type
       */
      template <typename Analysis>
      struct EdgeMaker :
         Loki::BaseVisitor,
         Loki::Visitor<GridPatch<CartesianSystem, Analysis>,
                       GridEdgePointer<Analysis>, true>,
         Loki::Visitor<GridPatch<PolarSystem, Analysis>,
                       GridEdgePointer<Analysis>, true>
      {
      private:
         using GridPatchCartesian =
            GridPatch<CartesianSystem, Analysis>;
         using GridPatchPolar =
            GridPatch<PolarSystem, Analysis>;
         using GridEdgePtr = GridEdgePointer<Analysis>;

      public:
         //! Initialization
         /*!
          *  \param row_index is the edge's starting vertex' row index.
          *
          *  \param column_index is the edge's starting vertex' column
          *  index.
          *
          *  \param direction is the edge's direction
          */
         EdgeMaker(
            RowIndex row_index, ColumnIndex column_index, Axis direction)
            : row_index_(row_index), column_index_(column_index),
              direction_(direction)
         {
         }

         EdgeMaker(const EdgeMaker&) = delete;

         GridEdgePtr Visit(const GridPatchCartesian&) override;
         GridEdgePtr Visit(const GridPatchPolar&) override;

      private:
         GridEdgePtr MakeStraightEdge(const GridPatchBase<Analysis>&) const;

      private:
         RowIndex row_index_;
         ColumnIndex column_index_;
         Axis direction_;
      };

      //--------------------------------------------------------------
      //
      //   Function implementations
      //
      //--------------------------------------------------------------

      template <typename Analysis>
      auto EdgeMaker<Analysis>::MakeStraightEdge(
         const GridPatchBase<Analysis>& grid) const -> GridEdgePtr
      {
         const auto edge_length = grid.EdgeLength(
            row_index_, column_index_, direction_);
         const auto crossing_edge_length = grid.CrossingEdgeLength(
            row_index_, column_index_, direction_);
         const auto reluctivity = grid.GetMaterial(
            row_index_, column_index_, direction_)->GetReluctivityModel(
               Flip(direction_));

         return std::make_unique<StraightGridEdge<Analysis> >(
            edge_length, crossing_edge_length/edge_length, reluctivity);
      }

      template <typename Analysis>
      auto EdgeMaker<Analysis>::Visit(
         const GridPatchCartesian& grid) -> GridEdgePtr
      {
         return this->MakeStraightEdge(grid);
      }

      template <typename Analysis>
      auto EdgeMaker<Analysis>::Visit(
         const GridPatchPolar& grid) -> GridEdgePtr
      {
         using std::get;

         if (Axis::kFirst == direction_)
            return this->MakeStraightEdge(grid);

         const auto edge_length = grid.EdgeLength(
            row_index_, column_index_, direction_);
         const auto crossing_edge_length = grid.CrossingEdgeLength(
            row_index_, column_index_, direction_);
         const Length radius = get<0>(
            grid.VertexPositionLocal(row_index_, column_index_));
         const auto reluctivity = grid.GetMaterial(
            row_index_, column_index_, direction_)->GetReluctivityModel(
               Flip(direction_));

         return std::make_unique<CircularGridEdge<Analysis> >(
            edge_length, crossing_edge_length/edge_length, radius,
            reluctivity);
      }
   }
}

#endif // EDGE_MAKER_HPP_Dengidat_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
