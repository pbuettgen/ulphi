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
 *  \brief Grid patch base class (geometry independent)
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef GRID_PATCH_BASE_HPP_STMC1ryi_
#define GRID_PATCH_BASE_HPP_STMC1ryi_

#include <set>

#include "config.h"

#include <boost/format.hpp>

#include "boundary.hpp"
#include "boundary_condition-fwd.hpp"
#include "coordinate_systems.hpp"
#include "grid.hpp"
#include "indices.hpp"
#include "load_definition-fwd.hpp"
#include "material-fwd.hpp"
#include "vertex_distribution.hpp"
#include "types.hpp"

namespace PACKAGE_NS {
   //! Geometry independent aspects of a grid patch
   /*!
    * @tparam Analysis is the analysis type
    */
   template<typename Analysis>
   struct GridPatchBase: Grid {
   protected:
      //! Base class
      using Base = Grid;

      //! A vertex's row and column index
      using RowColumnTuple = std::tuple<RowIndex, ColumnIndex>;

      //! Vertex index set
      using VertexIndexSetPointer = std::shared_ptr<std::set<VertexIndex> >;

      using BoundaryCondPointer = ConstBoundaryConditionPointer<Analysis>;
      using LoadDefPointer = ConstLoadDefinitionPointer<Analysis>;
      // using BoundaryAndVertex = std::tuple<const Boundary&, BoundaryIndex>;

   public:
      //! Destruction
      virtual ~GridPatchBase();

   protected:
      //! Construct a GridRegionBase
      /*!
       * @param vertex_distribution_horizontal is the vertex distribution for
       * the horizontal direction
       *
       * @param vertex_distribution_vertical is the vertex distribution for
       * the vertical direction
       */
      GridPatchBase(ConstVertexDistributionPointer,
                    ConstVertexDistributionPointer);

      //! Copy construction (unsupported)
      GridPatchBase(const GridPatchBase&) = delete;

   public:
      //! Assignment is unsupported
      GridPatchBase& operator = (const GridPatchBase&) = delete;

      //! Set a new vertex distribution for the specified axis
      /*!
       * @param vertex_distribution is the new vertex distribution
       *
       * @param axis is the axis for which the distribution shall be set
       */
      void SetVertexDistribution(ConstVertexDistributionPointer, Axis);

      //! Compute a vertex number from a vertex's row and column index
      /*!
       * Compute a vertex's global vertex number.  (Precondition is, that
       * SetFirstVertexNumber() has been called before.)  All indices
       * start counting with zero.
       *
       * The given row and column index is validated before a vertex
       * number is computed.  For boundary vertices this function takes
       * also into account that the boundary may have been coupled to
       * an other grid region's boundary.
       *
       * @param row_index is the vertex's row index
       * @param column_index is the vertex's column index
       * @return The computed vertex number
       * @throw OutOfRange in the case of an invalid row or column index.
       */
      VertexIndex GlobalIndex(
         RowIndex row_index, ColumnIndex column_index,
         VertexIndexSetPointer previous_vertices = nullptr) const;

      //! Compute a vertex number from a row-column-tuple
      /*!
       * Compute a vertex number without taking boundary couplings
       * into account.  For boundary vertices no forwarding to the
       * boundary is done.
       *
       * @param row_index is the vertex' row index
       * @param column_index is the vertex' column index
       * @return the vertex' index
       */
      VertexIndex GlobalIndexHere(
         RowIndex row_index, ColumnIndex column_index) const;

      //! Retrieve the number of vertex columns
      /*!
       * @return the number of vertex columns
       */
      VertexIndex NumberOfVertexColumns() const;

      //! Maximum vertex column index
      ColumnIndex MaxColumnIndex() const {
         return ColumnIndex(NumberOfVertexColumns()-1);
      }

      //! Retrieve the number of vertex rows
      /*!
       * @return the number of vertex rows.
       */
      VertexIndex NumberOfVertexRows() const;

      //! Maximum vertex row index
      RowIndex MaxRowIndex() const {
         return RowIndex(NumberOfVertexRows()-1);
      }

      //! Compute the number of vertices within this grid region
      /*!
       * Also count vertices on coupled boundaries.  This number is
       * relevant to predict the maximum vertex index.
       *
       * @return the number of vertices within this grid patch.
       */
      VertexIndex NumberOfVerticesBrutto() const override;

      //! Compute the number of vertices within this grid
      /*!
       * Ignore vertices on coupled  boundaries.
       *
       * @return total number of vertices - number of vertices on
       * coupled boundaries.
       */
      VertexIndex NumberOfVerticesNetto() const;

      //! Compute the number of edges within this grid
      /*!
       * Also count edges on coupled boundaries
       *
       * @return number of edges within this grid
       */
      VertexIndex NumberOfEdgesBrutto() const;

      //! Compute the number of edges within this grid
      /*!
       * Substracts the number of edges on coupled boundaries from
       * total number of edges.
       *
       * @return number of edges (brutto) - number of edges on coupled
       * boundaries
       */
      VertexIndex NumberOfEdgesNetto() const override;

      //! Decompose a vertex index into its RowColumnTuple
      /*!
       * @param i is the VertexIndex to decompose. It is considered as
       * a local vertex index.
       *
       * @return a RowColumnTuple corresponding to the vertex's row
       * and column index.
       */
      RowColumnTuple RowColumnDecompose(VertexIndex) const;

      //! Compute the edge length
      /*!
       * @param row_index the source vertex' row index
       *
       * @param column_index the source vertex' column index
       *
       * @param axis direction
       */
      virtual Length EdgeLength(RowIndex, ColumnIndex, Axis) const =0;

      //! Average edge length at a given vertex
      /*!
       * @param row_index is the vertex' row index
       *
       * @param column_index is the vertex' column index
       *
       * @param axis direction
       */
      Length AverageEdgeLength(RowIndex, ColumnIndex, Axis) const;

      //! Get load at a specified vertex
      /*!
       * @param row_index is the vertex' row index
       *
       * @param column_index is the vertex' column index
       *
       * @returns the load condition.
       */
      LoadDefPointer GetLoad(RowIndex, ColumnIndex) const;

      //! Get the material at a specified vertex
      /*!
       * @param row_index is the vertex' row index
       *
       * @param column_index is the vertex' column index
       *
       * @return a pointer to the material if a material for the
       * specified point is found or nullpointer else.
       */
      ConstMaterialPointer GetMaterial(RowIndex, ColumnIndex) const;

      //! Get the material at a edge'es midpoint
      /*!
       * The edge is described by its starting vertex and its direction.
       *
       * @param row_index is the edge'es starting vertex' row index.
       *
       * @param column_index is the edge'es starting vertex' column
       * index.
       *
       * @return a pointer to the material if a material for the
       * specified edge is found or nullpointer else.
       */
      ConstMaterialPointer GetMaterial(RowIndex, ColumnIndex, Axis) const;

      //! Get the material at a specified vertex
      /*!
       * Only look at this grid; no forwarding to the boundary.
       *
       * @param row_index is the vertex' row index.
       *
       * @param column_index is the vertex' column index.
       *
       * @returns a pointer to the material at the specified vertex.
       */
      virtual ConstMaterialPointer GetMaterialHere(
         RowIndex, ColumnIndex) const =0;

      //! Get the material at a specified edge
      /*!
       * Only look at this grid; no forwarding to the boundary. The
       * material is determined at the edge'es midpoint.
       *
       * @param row_index specifies the edge'es source vertex
       *
       * @param column_index specifies the edge'es source vertex
       *
       * @param axis is the edge'es direction.
       *
       * @returns a pointer to the material at the specified vertex.
       */
      virtual ConstMaterialPointer GetMaterialHere(
         RowIndex, ColumnIndex, Axis) const =0;

      //! Get the load definition at a specified vertex
      /*!
       * For boundary vertices no forwarding to the boundary is done.
       *
       * @param row_index the vertex' row index
       *
       * @param column_index the vertex' column index
       */
      virtual LoadDefPointer GetLoadHere(RowIndex, ColumnIndex) const =0;

      //! Length of the secondary grid's crossing edge
      /*!
       * The parameters specify the primary grid's edge.
       *
       * @param row_index the vertex' row index
       *
       * @param column_index the vertex' column index
       *
       * @param axis is the edge'es direction.
       */
      Length CrossingEdgeLength(RowIndex, ColumnIndex, Axis) const;

      //! Length of the secondary grid's crossing edge
      /*!
       * For boundary vertices no forwarding to the boundary is done.
       *
       * @param row_index the vertex' row index
       *
       * @param column_index the vertex' column index
       *
       * @param axis is the edge'es direction.
       */
      virtual Length CrossingEdgeLengthHere(
         RowIndex, ColumnIndex, Axis) const =0;

      //! Secondary grid's cell size surround a vertex
      /*!
       * @param row_index the vertex' row index
       *
       * @param column_index the vertex' column index
       */
      Area CellSize(RowIndex, ColumnIndex) const;

      //! Secondary grid's cell size surround a vertex
      /*!
       * For boundary vertices no forwarding to the boundary is done.
       *
       * @param row_index the vertex' row index
       *
       * @param column_index the vertex' column index
       */
      virtual Area CellSizeHere(RowIndex, ColumnIndex) const =0;

      //! Compute the bottom half of a cell surround a vertex
      virtual Area CellSizeBottomHalf(RowIndex, ColumnIndex) const =0;

      //! Compute the top half of a cell surround a vertex
      virtual Area CellSizeTopHalf(RowIndex, ColumnIndex) const =0;

      //! Compute the left half of a cell surround a vertex
      virtual Area CellSizeLeftHalf(RowIndex, ColumnIndex) const =0;

      //! Compute the right half of a cell surround a vertex
      virtual Area CellSizeRightHalf(RowIndex, ColumnIndex) const =0;

      virtual CartesianSystem::PositionVector VertexPositionGlobal(
         RowIndex, ColumnIndex) const =0;

      //! Check that all coupled boundaries are compatible
      void CheckBoundaryCompability() const;

      //! Get a reference to the grid's east edge
      /*!
       * @return a reference to the grid's east edge
       */
      EastBoundary<Analysis>& east_boundary();

      //! Get a reference to the grid's east edge
      /*!
       * @return a reference to the grid's east edge
       */
      const EastBoundary<Analysis>& east_boundary() const;

      //! Get a reference to the grid's north edge
      /*!
       * @return a reference to the grid's north edge
       */
      NorthBoundary<Analysis>& north_boundary();

      //! Get a reference to the grid's north edge
      /*!
       * @return a reference to the grid's north edge
       */
      const NorthBoundary<Analysis>& north_boundary() const;

      //! Get a reference to the grid's west edge
      /*!
       * @return a reference to the grid's west edge
       */
      WestBoundary<Analysis>& west_boundary();

      //! Get a reference to the grid's west edge
      /*!
       * @return a reference to the grid's west edge
       */
      const WestBoundary<Analysis>& west_boundary() const;

      //! Get a reference to the grid's south edge
      /*!
       * @return a reference to the grid's south edge
       */
      SouthBoundary<Analysis>& south_boundary();

      //! Get a reference to the grid's south edge
      /*!
       * @return a reference to the grid's south edge
       */
      const SouthBoundary<Analysis>& south_boundary() const;

   protected:
      //! Hook for assigning the first vertex index in the global index space
      /*!
       * @param vertex_index base index used by this grid
       * @return next free vertex index in the global index space
       */
      VertexIndex DoSetFirstVertexIndex(
         VertexIndex vertex_index) const override;

      //! Check whether the given row index is valid
      /*!
       * @param row is the row index to check
       * @throw OutOfRange in the case of an invalid vertex row index
       */
      void CheckValid(RowIndex) const;

      //! Check whether the give column index is valid
      /*!
       * @param column is the column index to check
       * @throw OutOfRange in the case of an invalid vertex column index
       */
      void CheckValid(ColumnIndex) const;

      void CheckNoCoupling(const Boundary<Analysis>&) const;

      //! Compute a vertex's normalized position within a row
      /*!
       * @param column is the vertex's column index
       * @return the vertex's normalized position within a row
       */
      double VertexPositionInRow(ColumnIndex) const;

      //! Compute a vertex's normalized position within a column
      /*!
       * @param row is the vertex's row index
       * @return the vertex's normalized position within a column
       */
      double VertexPositionInColumn(RowIndex) const;

   private:
      bool CheckCyclicBoundaryIndexComputation(VertexIndex,
                                               VertexIndexSetPointer&) const;

   private:
      EastBoundary<Analysis> east_boundary_;
      NorthBoundary<Analysis> north_boundary_;
      WestBoundary<Analysis> west_boundary_;
      SouthBoundary<Analysis> south_boundary_;
      ConstVertexDistributionPointer vertex_distribution_horizontal_,
         vertex_distribution_vertical_;
   };

   //--------------------------------------------------------------
   //
   //   Function implementations
   //
   //--------------------------------------------------------------

   template<typename Analysis>
   inline GridPatchBase<Analysis>::GridPatchBase(
      ConstVertexDistributionPointer vertex_distribution_horizontal,
      ConstVertexDistributionPointer vertex_distribution_vertical)
      : east_boundary_(this), north_boundary_(this),
        west_boundary_(this), south_boundary_(this),
        vertex_distribution_horizontal_(vertex_distribution_horizontal),
        vertex_distribution_vertical_(vertex_distribution_vertical)
   {
   }

   template<typename Analysis>
   GridPatchBase<Analysis>::~GridPatchBase()
   {
   }

   template <typename Analysis>
   inline void GridPatchBase<Analysis>::SetVertexDistribution(
      ConstVertexDistributionPointer vertex_distribution, Axis axis)
   {
      if (Axis::kFirst == axis)
         vertex_distribution_horizontal_ = vertex_distribution;
      else
         vertex_distribution_vertical_ = vertex_distribution;
   }

   template<typename Analysis>
   inline VertexIndex GridPatchBase<Analysis>::DoSetFirstVertexIndex(
      VertexIndex vertex_number) const
   {
      return this->ConvertLocal2Global(NumberOfVerticesBrutto());
   }

   template <typename Analysis>
   inline Length GridPatchBase<Analysis>::CrossingEdgeLength(
      RowIndex row_index, ColumnIndex column_index, Axis axis) const
   {
#ifndef NDEBUG
      CheckValid(row_index);
      CheckValid(column_index);
      if (Axis::kFirst == axis)
         this->CheckValid(column_index + 1);
      else
         this->CheckValid(row_index + 1);
#endif

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      if (Axis::kFirst == axis) {
         if (0 == row_index) {
            return south_boundary_.CrossingEdgeLength(column_index);
         }
         if (max_row_index == row_index) {
            return north_boundary_.CrossingEdgeLength(column_index);
         }
      }
      else {
         if (0 == column_index) {
            return west_boundary_.CrossingEdgeLength(row_index);
         }
         if (max_column_index == column_index) {
            return east_boundary_.CrossingEdgeLength(row_index);
         }
      }

      return CrossingEdgeLengthHere(row_index, column_index, axis);
   }

   template <typename Analysis>
   inline Length GridPatchBase<Analysis>::AverageEdgeLength(
      RowIndex row_index, ColumnIndex column_index, Axis axis) const
   {
#ifndef NDEBUG
      this->CheckValid(row_index);
      this->CheckValid(column_index);
#endif

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      if (Axis::kFirst == axis) {
         if (0 == column_index)
            return west_boundary().AverageEdgeLengthPerpendicular(row_index);
         if (max_column_index == column_index)
            return east_boundary().AverageEdgeLengthPerpendicular(row_index);
      }
      else {
         if (0 == row_index)
            return south_boundary().AverageEdgeLengthPerpendicular(
               column_index);
         if (max_row_index == row_index)
            return north_boundary().AverageEdgeLengthPerpendicular(
               column_index);
      }

      return (this->EdgeLength(row_index, column_index, axis)
              + this->EdgeLength(
                 row_index - (Axis::kSecond == axis),
                 column_index - (Axis::kFirst == axis), axis))/2.;
   }

   template <typename Analysis>
   inline auto GridPatchBase<Analysis>::GetLoad(
      RowIndex row_index, ColumnIndex column_index) const -> LoadDefPointer
   {
#ifndef NDEBUG
      this->CheckValid(row_index);
      this->CheckValid(column_index);
#endif

      using std::get;

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      if (0==row_index) {
         if (0==column_index) {
            // South west corner
            return this->south_boundary().HasBoundaryCoupling() ?
               this->south_boundary().GetLoad(column_index):
               this->west_boundary().GetLoad(row_index);
         }
         if (max_column_index == column_index) {
            // South east corner
            return this->south_boundary().HasBoundaryCoupling() ?
               this->south_boundary().GetLoad(column_index):
               this->east_boundary().GetLoad(row_index);
         }
         // South boundary
         return this->south_boundary().GetLoad(column_index);
      }
      if (max_row_index == row_index) {
         if (0==column_index) {
            // North west corner
            return this->north_boundary().HasBoundaryCoupling() ?
               this->north_boundary().GetLoad(column_index):
               this->west_boundary().GetLoad(row_index);
         }
         if (max_column_index == column_index) {
            // North east corner
            return this->north_boundary().HasBoundaryCoupling() ?
               this->north_boundary().GetLoad(column_index):
               this->east_boundary().GetLoad(row_index);
         }
         // North boundary
         return this->north_boundary().GetLoad(column_index);
      }
      if (0 == column_index) {
         // West boundary
         return this->west_boundary().GetLoad(row_index);
      }
      if (max_column_index == column_index) {
         // East boundary
         return this->east_boundary().GetLoad(row_index);
      }

      return this->GetLoadHere(row_index, column_index);
   }

   template <typename Analysis>
   inline ConstMaterialPointer GridPatchBase<Analysis>::GetMaterial(
      RowIndex row_index, ColumnIndex column_index, Axis axis) const
   {
#ifndef NDEBUG
      CheckValid(row_index);
      CheckValid(column_index);
#endif

      using std::get;

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      if (Axis::kFirst == axis) {
         if (0 == row_index) {
            // South boundary
            return south_boundary().GetMaterialAtEdge(column_index);
         }
         if (max_row_index == row_index) {
            // North boundary
            return north_boundary().GetMaterialAtEdge(column_index);
         }
      }
      else {
         if (0 == column_index) {
            // West boundary
            return west_boundary().GetMaterialAtEdge(row_index);
         }
         if (max_column_index == column_index) {
            // East boundary
            return east_boundary().GetMaterialAtEdge(row_index);
         }
      }

      return GetMaterialHere(row_index, column_index, axis);
   }

   template <typename Analysis>
   inline ConstMaterialPointer GridPatchBase<Analysis>::GetMaterial(
      RowIndex row_index, ColumnIndex column_index) const
   {
#ifndef NDEBUG
      CheckValid(row_index);
      CheckValid(column_index);
#endif

      using std::get;

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      if (0 == row_index) {
         if (0 == column_index) {
            // South west corner
            return this->SouthHasBoundaryCoupling() ?
               south_boundary().GetMaterialAtVertex(column_index):
               west_boundary().GetMaterialAtVertex(row_index);
         }
         if (max_column_index == column_index) {
            // South east corner
            return this->SouthHasBoundaryCoupling() ?
               south_boundary().GetMaterialAtVertex(column_index) :
               east_boundary().GetMaterialAtVertex(row_index);
         }
         // South boundary
         return south_boundary().GetMaterialAtVertex(column_index);
      }
      if (max_row_index == row_index) {
         if (0==column_index) {
            // North west corner
            return this->NorthHasBoundaryCoupling() ?
               north_boundary().GetMaterialAtVertex(column_index) :
               west_boundary().GetMaterialAtVertex(row_index);
         }
         if (max_column_index == column_index) {
            // North east corner
            return this->NorthHasBoundaryCoupling() ?
               north_boundary().GetMaterialAtVertex(column_index) :
               east_boundary().GetMaterialAtVertex(row_index);
         }
         // North boundary
         return north_boundary().GetMaterialAtVertex(column_index);
      }
      if (0 == column_index) {
         // West boundary
         return west_boundary().GetMaterialAtVertex(row_index);
      }
      if (max_column_index == column_index) {
         return east_boundary().GetMaterialAtVertex(row_index);
      }

      return this->GetMaterialHere(row_index, column_index);
   }

   template<typename Analysis>
   inline VertexIndex GridPatchBase<Analysis>::NumberOfVertexColumns() const
   {
      return vertex_distribution_horizontal_->number_of_vertices();
   }

   template<typename Analysis>
   inline VertexIndex GridPatchBase<Analysis>::NumberOfVertexRows() const
   {
      return vertex_distribution_vertical_->number_of_vertices();
   }

   template<typename Analysis>
   inline VertexIndex GridPatchBase<Analysis>::NumberOfVerticesNetto() const
   {
      VertexIndex number_of_vertices = this->NumberOfVerticesBrutto();

      if (north_boundary_.HasBoundaryCoupling())
         number_of_vertices -= north_boundary_.NumberOfVertices();
      if (west_boundary_.HasBoundaryCoupling())
         number_of_vertices -= west_boundary_.NumberOfVertices();
      if (south_boundary_.HasBoundaryCoupling())
         number_of_vertices -= south_boundary_.NumberOfVertices();
      if (east_boundary_.HasBoundaryCoupling())
         number_of_vertices -= east_boundary_.NumberOfVertices();

      return number_of_vertices;
   }

   template<typename Analysis>
   inline VertexIndex GridPatchBase<Analysis>::NumberOfVerticesBrutto() const
   {
      return NumberOfVertexRows() * NumberOfVertexColumns();
   }

   template<typename Analysis>
   inline VertexIndex GridPatchBase<Analysis>::NumberOfEdgesBrutto() const
   {
      const VertexIndex number_of_vertex_rows = this->NumberOfVertexRows();
      const VertexIndex number_of_vertex_colums = this->NumberOfVertexColumns();

      return 2 * number_of_vertex_rows * number_of_vertex_colums
         - number_of_vertex_rows - number_of_vertex_colums;
   }

   template<typename Analysis>
   inline VertexIndex GridPatchBase<Analysis>::NumberOfEdgesNetto() const
   {
      VertexIndex number_of_edges = this->NumberOfEdgesBrutto();

      if (north_boundary_.HasBoundaryCoupling())
         number_of_edges -= north_boundary_.NumberOfEdges();
      if (west_boundary_.HasBoundaryCoupling())
         number_of_edges -= west_boundary_.NumberOfEdges();
      if (south_boundary_.HasBoundaryCoupling())
         number_of_edges -= south_boundary_.NumberOfEdges();
      if (east_boundary_.HasBoundaryCoupling())
         number_of_edges -= east_boundary_.NumberOfEdges();

      return number_of_edges;
   }

   template <typename Analysis>
   inline void GridPatchBase<Analysis>::CheckNoCoupling(
      const Boundary<Analysis>& boundary) const {
      if (boundary.HasBoundaryCoupling()) {
         boost::format error_message(
            "%1%: Attemp to create invalid coupling pattern!");
#ifdef __GNUC__
#define FUNCTION_NAME __PRETTY_FUNCTION__
#else
#define FUNCTION_NAME typeid(*this).name()
#endif
         throw OutOfRange((error_message % FUNCTION_NAME).str());
#undef FUNCTION_NAME
      }
   }

   template<typename Analysis>
   inline void GridPatchBase<Analysis>::CheckValid(
      RowIndex row_index) const
   {
      const auto max_row_index = this->MaxRowIndex();

      if (max_row_index < row_index) {
         boost::format error_message("%1%: Invalid row index %2% (max: %3%)!");
#ifdef __GNUC__
#define FUNCTION_NAME __PRETTY_FUNCTION__
#else
#define FUNCTION_NAME typeid(*this).name()
#endif
         throw OutOfRange(
            (error_message % FUNCTION_NAME % row_index % max_row_index).str());
#undef FUNCTION_NAME
      }
   }

   template<typename Analysis>
   inline void GridPatchBase<Analysis>::CheckValid(
      ColumnIndex column_index) const
   {
      const auto max_column_index = this->MaxColumnIndex();

      if (max_column_index < column_index) {
         boost::format error_message(
            "%1%: Invalid column index %2% (max: %3%)!");

#ifdef __GNUC__
#define FUNCTION_NAME __PRETTY_FUNCTION__
#else
#define FUNCTION_NAME typeid(*this).name()
#endif
         throw OutOfRange(
            (error_message % FUNCTION_NAME % column_index
             % max_column_index).str());
#undef FUNCTION_NAME
      }
   }

   template<typename Analysis>
   inline void GridPatchBase<Analysis>::CheckBoundaryCompability() const
   {
      east_boundary_.CheckYourPartner();
      north_boundary_.CheckYourPartner();
      west_boundary_.CheckYourPartner();
      south_boundary_.CheckYourPartner();
   }

   template<typename Analysis>
   bool GridPatchBase<Analysis>::CheckCyclicBoundaryIndexComputation(
      VertexIndex global_vertex_index,
      VertexIndexSetPointer& previous_vertices) const
   {
      if (nullptr != previous_vertices) {
         if (previous_vertices->count(global_vertex_index))
            return true;
      } else {
         previous_vertices = std::make_shared<std::set<VertexIndex> >();
      }
      previous_vertices->insert(global_vertex_index);
      return false;
   }

   template<typename Analysis>
   VertexIndex GridPatchBase<Analysis>::GlobalIndex(
      RowIndex row_index, ColumnIndex column_index,
      VertexIndexSetPointer previous_vertices) const
   {
#ifndef NDEBUG
      CheckValid(row_index);
      CheckValid(column_index);
#endif

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      if (0 == row_index)
      {
         if (0 == column_index) {
            // south west corner
            const auto global_vertex_index
               = this->GlobalIndexHere(row_index, column_index);
            if (CheckCyclicBoundaryIndexComputation(
                   global_vertex_index, previous_vertices))
               return global_vertex_index;

            return west_boundary_.HasBoundaryCoupling() ?
               west_boundary_.GlobalIndex(row_index, previous_vertices) :
               south_boundary_.GlobalIndex(column_index, previous_vertices);
         }
         if (max_column_index == column_index)
         {
            // south east corner
            const auto global_vertex_index
               = this->GlobalIndexHere(row_index, column_index);
            if (CheckCyclicBoundaryIndexComputation(
                   global_vertex_index, previous_vertices))
               return global_vertex_index;

            return south_boundary_.HasBoundaryCoupling() ?
               south_boundary_.GlobalIndex(column_index, previous_vertices) :
               east_boundary_.GlobalIndex(row_index, previous_vertices);
         }
         // south boundary
         return south_boundary_.GlobalIndex(column_index, previous_vertices);
      }
      if (max_row_index == row_index)
      {
         if (max_column_index == column_index) {
            // north east corner
            const auto global_vertex_index
               = this->GlobalIndexHere(row_index, column_index);
            if (CheckCyclicBoundaryIndexComputation(
                   global_vertex_index, previous_vertices))
               return global_vertex_index;

            return east_boundary_.HasBoundaryCoupling() ?
               east_boundary_.GlobalIndex(row_index, previous_vertices) :
               north_boundary_.GlobalIndex(column_index, previous_vertices);
         }
         if (0 == column_index)
         {
            // north west corner
            const auto global_vertex_index
               = this->GlobalIndexHere(row_index, column_index);
            if (CheckCyclicBoundaryIndexComputation(
                   global_vertex_index, previous_vertices))
               return global_vertex_index;

            return north_boundary_.HasBoundaryCoupling() ?
               north_boundary_.GlobalIndex(column_index, previous_vertices):
               west_boundary_.GlobalIndex(row_index, previous_vertices);
         }
         // north boundary
         return north_boundary_.GlobalIndex(column_index, previous_vertices);
      }
      if (max_column_index == column_index)
      {
         // east boundary
         return east_boundary_.GlobalIndex(row_index, previous_vertices);
      }
      if (0 == column_index)
      {
         // west boundary
         return west_boundary_.GlobalIndex(row_index, previous_vertices);
      }

      return GlobalIndexHere(row_index, column_index);
   }

   template<typename Analysis>
   inline VertexIndex GridPatchBase<Analysis>::GlobalIndexHere(
      RowIndex row_index, ColumnIndex column_index) const
   {
      return this->ConvertLocal2Global(
         column_index.Get() + NumberOfVertexColumns() * row_index.Get());
   }

   template <typename Analysis>
   inline Area GridPatchBase<Analysis>::CellSize(
      RowIndex row_index, ColumnIndex column_index) const
   {
#ifndef NDEBUG
      this->CheckValid(column_index);
      this->CheckValid(row_index);
#endif
      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      if (0 == row_index) {
         if (0 == column_index) {
            // South west corner
            return this->south_boundary().HasBoundaryCoupling() ?
               this->south_boundary().CellSize(column_index):
               this->west_boundary().CellSize(row_index);
         }
         if (max_column_index == column_index) {
            // South east corner
            return this->south_boundary().HasBoundaryCoupling() ?
               this->south_boundary().CellSize(column_index):
               this->east_boundary().CellSize(row_index);
         }
         // South boundary
         return this->south_boundary().CellSize(column_index);
      }
      if (max_row_index == row_index) {
         if (0 == column_index) {
            // North west corner
            return this->north_boundary().HasBoundaryCoupling() ?
               this->north_boundary().CellSize(column_index):
               this->west_boundary().CellSize(row_index);
         }
         if (max_column_index == column_index) {
            // North east corner
            return this->north_boundary().HasBoundaryCoupling() ?
               this->north_boundary().CellSize(column_index):
               this->east_boundary().CellSize(row_index);
         }
         // North boundary
         return this->north_boundary().CellSize(column_index);
      }
      if (0 == column_index) {
         // West boundary
         return this->west_boundary().CellSize(row_index);
      }
      if (max_column_index == column_index) {
         // East boundary
         return this->east_boundary().CellSize(row_index);
      }

      return CellSizeHere(row_index, column_index);
   }

   template<typename Analysis>
   inline auto GridPatchBase<Analysis>::RowColumnDecompose(
      VertexIndex local_vertex_index) const -> RowColumnTuple
   {
      const auto number_of_vertex_columns = NumberOfVertexColumns();
      return std::make_tuple(
         RowIndex(local_vertex_index / number_of_vertex_columns),
         ColumnIndex(local_vertex_index % number_of_vertex_columns));
   }

   template<typename Analysis>
   inline double GridPatchBase<Analysis>::VertexPositionInRow(
      ColumnIndex column_index) const
   {
      return vertex_distribution_horizontal_->PositionNthVertex(
         column_index.Get());
   }

   template<typename Analysis>
   inline double GridPatchBase<Analysis>::VertexPositionInColumn(
      RowIndex row_index) const
   {
      return vertex_distribution_vertical_->PositionNthVertex(
         row_index.Get());
   }

#define ACCESS_BOUNDARY(DIR1, DIR2)                                     \
   template<typename Analysis>                                          \
   inline DIR1 ## Boundary<Analysis>& GridPatchBase<                    \
                                                        Analysis>:: DIR2 ## _boundary() \
   {                                                                    \
      return DIR2 ## _boundary_;                                        \
   }                                                                    \
   template<typename Analysis>                                          \
   inline const DIR1 ## Boundary<Analysis>& GridPatchBase<              \
                                                                        Analysis>:: DIR2 ## _boundary() const \
   {                                                                    \
      return DIR2 ## _boundary_;                                        \
   }

   ACCESS_BOUNDARY(East, east)
   ACCESS_BOUNDARY(North, north)
   ACCESS_BOUNDARY(West, west)
   ACCESS_BOUNDARY(South, south)
#undef ACCESS_BOUNDARY
}

#endif // GRID_PATCH_BASE_HPP_STMC1ryi_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
