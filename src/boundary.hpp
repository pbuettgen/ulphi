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
 *  \brief A GridPatch's boundary
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef BOUNDARY_HPP_FGr2WLZN_
#define BOUNDARY_HPP_FGr2WLZN_

#include "config.h"

#include <boost/format.hpp>

#include "analysis.hpp"
#include "coordinate_systems.hpp"
#include "boundary_condition-fwd.hpp"
#include "exception.hpp"
#include "grids-fwd.hpp"
#include "indices.hpp"
#include "load_definition-fwd.hpp"
#include "material-fwd.hpp"
#include "types.hpp"
#include "void_condition.hpp"

namespace PACKAGE_NS {
   //! Represent a GridPatch's boundary
   /*!
    * @tparam Analysis is the selected analysis type
    */
   template<typename Analysis>
   struct Boundary {
   protected:
      //! Pointer to an boundary object
      using ConstBoundaryPointer = const Boundary<Analysis>*;

      //! Pointer to the grid this boundary belongs to
      using ConstGridPatchBasePointer = const GridPatchBase<Analysis>*;

      //! A set of vertex indices
      using VertexIndexSetPointer = std::shared_ptr<std::set<VertexIndex> >;

      //! A vertex's row and column index
      using VertexID = std::tuple<RowIndex, ColumnIndex>;

      //! An edge specifier
      using EdgeID = std::tuple<RowIndex, ColumnIndex, Axis>;

      //! A pointer to a LoadDefinition object
      using LoadDefPointer = ConstLoadDefinitionPointer<Analysis>;

   protected:
      //! Construct a boundary
      /*!
       * @param grid_patch points to the GridPatch this edge belongs
       * to.
       */
      explicit Boundary(ConstGridPatchBasePointer grid_patch) :
         grid_patch_(grid_patch),
         boundary_condition_(nullptr),
         default_condition_(VoidCondition<Analysis>::New()),
         partner_boundary_(nullptr)
      {
      }

   public:
      //! Destruct a boundary
      virtual ~Boundary()
      {
      }

      //! Assignment is unsupported
      Boundary& operator =(const Boundary&) = delete;

      //! Compute the number of vertices onto this boundary
      virtual VertexIndex NumberOfVertices() const =0;

      //! Compute the number of grid edges on this boundary
      virtual VertexIndex NumberOfEdges() const =0;

      //! Reverse a vertex index
      /*!
       * @param i is the vertex index to reverse
       *
       * @return the reversed vertex index
       */
      BoundaryIndex ReverseIndex(BoundaryIndex i) const;

      //! Return the global index for a given vertex onto the boundary
      /*!
       * All indices start with zero.
       *
       * @param i is the vertex' index onto the boundary.
       *
       * @param previous_vertices is used to detect cycles in the
       * vertex number computation.  (This may happen for vertices at
       * the corners.)
       *
       * @return the vertex' global index.
       *
       * @throw OutOfRange in the case of an invalid index.
       */
      VertexIndex GlobalIndex(BoundaryIndex, VertexIndexSetPointer) const;

      //! Get the corresponding row and column index.
      /*!
       * Get the corresponding row and column index for a given index
       * onto the boundary.
       *
       * @returns a row and column tuple
       */
      virtual VertexID AsGridVertexID(BoundaryIndex) const =0;

      //! Compute a EdgeID from a BoundaryIndex
      virtual EdgeID AsGridEdgeID(BoundaryIndex) const =0;

      //! Grid spacing perpendicular to the boundary
      virtual Length EdgeLengthPerpendicular(BoundaryIndex) const =0;

      //! Average edge length
      Length AverageEdgeLengthPerpendicular(BoundaryIndex) const;

      //! Secondary grid's edge length
      Length CrossingEdgeLength(BoundaryIndex) const;

      //! Compute the secondary grid's cell size
      /*!
       * Compute the secondary grid's cell size suround a given vertex
       *
       * @param i is the vertex index
       *
       * @return cell size at given vertex index
       *
       * @bug The case of more than two grids sharing a corner is not
       * yet handled.
       */
      Area CellSize(BoundaryIndex) const;

      virtual Area CellSizeInnerHalf(BoundaryIndex) const =0;

      CartesianSystem::Point VertexPositionGlobal(BoundaryIndex) const;

      //! Check the partner edge for compability
      void CheckYourPartner() const;

      //! Couple this boundary to another boundary
      /*!
       * When this boundary is coupled to an other one, the other
       * boundary hides this boundary.
       *
       * @param other is the boundary which shall be coupled to this
       * boundary.
       */
      void SetBoundaryCoupling(const Boundary<Analysis>& other);

      //! Delete a previously defined boundary coupling
      void DeleteBoundaryCoupling();

      //! Determine whether this edge has a boundary coupling
      /*!
       * @return true if this edge is linked to an other edge
       */
      bool HasBoundaryCoupling() const;

      //! Set a boundary condition to this edge
      /*!
       * @param condition is the boundary condition to assign
       */
      void SetBoundaryCondition(ConstBoundaryConditionPointer<Analysis>);

      //! Delete a previously assigned boundary condition
      void DeleteBoundaryCondition();

      //! Determine whether this Boundary has a boundary condition set
      bool HasBoundaryCondition() const;

      //! Retrieve the boundary's condition
      /*!
       * @returns a previously set boundary condition or a VoidCondition
       */
      ConstBoundaryConditionPointer<Analysis> boundary_condition() const
      {
         return nullptr == boundary_condition_ ?
            default_condition_ : boundary_condition_;
      }

      //! Get the load condition at a specified vertex
      /*!
       * The operation is forwarded to the partner boundary (if any)
       * or routed back to the boundary's grid else.
       *
       * @param vertex is the vertex at which the load condition shall
       * be determined.
       *
       * @returns the load condition at vertex
       */
      LoadDefPointer GetLoad(BoundaryIndex) const;

      //! Get the material at a specified vertex
      /*!
       * The operation is forwarded to the partner boundary (if any)
       * or routed back to the boundary's grid else.
       *
       * @param vertex is the vertex at which the material shall
       * be determined.
       *
       * @returns the material at vertex
       */
      ConstMaterialPointer GetMaterialAtVertex(BoundaryIndex) const;

      //! Get the material at a specified edge
      /*!
       * The material at the edge'es midpoint is determined.
       *
       * The operation is forwarded to the partner boundary (if any)
       * or routed back to the boundary's grid else.
       *
       * @param vertex is the starting vertex of the edge for which
       * the material shall be determined.
       *
       * @returns the material at the edge
       */
      ConstMaterialPointer GetMaterialAtEdge(BoundaryIndex) const;

      virtual Parallelity parallelity() const =0;

   protected:
      //! Retrieve the grid to which this edge belongs
      /*!
       * @return the grid this edge belongs to
       */
      ConstGridPatchBasePointer grid_patch() const
      {
         return grid_patch_;
      }

      //! Access the boundary's partner boundary (if any)
      ConstBoundaryPointer partner_boundary() const
      {
         return partner_boundary_;
      }

   private:
      ConstGridPatchBasePointer grid_patch_;
      ConstBoundaryConditionPointer<Analysis> boundary_condition_;
      ConstBoundaryConditionPointer<Analysis> default_condition_;
      ConstBoundaryPointer partner_boundary_;
   };

   //! Function definitions useful for all boundaries
   /*!
    * These functions operate on a single vertex
    */
#define BOUNDARY_VERTEX_FUNCTIONS(INDEX_TYPE)                           \
   Length AverageEdgeLengthPerpendicular(INDEX_TYPE i) const {          \
      return Base::AverageEdgeLengthPerpendicular(AsBoundaryIndex(i));	\
   }                                                                    \
   VertexIndex GlobalIndex(INDEX_TYPE i, VertexIndexSetPointer p) const	\
   {                                                                    \
      return Base::GlobalIndex(AsBoundaryIndex(i), p);                  \
   }                                                                    \
   ConstMaterialPointer GetMaterialAtVertex(INDEX_TYPE i) const {       \
      return Base::GetMaterialAtVertex(AsBoundaryIndex(i));             \
   }                                                                    \
   LoadDefPointer GetLoad(INDEX_TYPE i) const {                         \
      return Base::GetLoad(AsBoundaryIndex(i));                         \
   }                                                                    \
   Area CellSize(INDEX_TYPE i) const {                                  \
      return Base::CellSize(AsBoundaryIndex(i));                        \
   }

   //! Function definitions usefull for all boundaries
   /*!
    * These functions operate on an edge.
    */
#define BOUNDARY_EDGE_FUNCTIONS(INDEX_TYPE, OFFSET)             \
   ConstMaterialPointer GetMaterialAtEdge(INDEX_TYPE i) const { \
   return Base::GetMaterialAtEdge(AsBoundaryIndex(i+OFFSET));	\
}                                                               \
   Length CrossingEdgeLength(INDEX_TYPE i) const {              \
   return Base::CrossingEdgeLength(AsBoundaryIndex(i+OFFSET));	\
}

   //! Declaration of virtual functions
#define BOUNDARY_VIRTUAL_FUNCTIONS                                  \
   VertexIndex NumberOfVertices() const override;                   \
   VertexIndex NumberOfEdges() const override;                      \
   Length EdgeLengthPerpendicular(BoundaryIndex) const override;	\
   VertexID AsGridVertexID(BoundaryIndex) const override;           \
   EdgeID AsGridEdgeID(BoundaryIndex) const override;               \
   Area CellSizeInnerHalf(BoundaryIndex) const override;            \
   Parallelity parallelity() const override;

   //! Types used by all boundaries
#define BOUNDARY_USING                                      \
   using Base = Boundary<Analysis>;                         \
   using ConstGridPatchBasePointer =                        \
      typename Base::ConstGridPatchBasePointer;             \
   using VertexIndexSetPointer =                            \
      typename Base::VertexIndexSetPointer;                 \
   using VertexID = typename Base::VertexID;                \
   using EdgeID = typename Base::EdgeID;                    \
   using LoadDefPointer = typename Base::LoadDefPointer;

   //! Represent a grid's north edge
   /*!
    * @tparam Analysis is the analysis type
    */
   template<typename Analysis>
   struct NorthBoundary: Boundary<Analysis> {
   private:
      BOUNDARY_USING

      public:
      //! Initialization
      /*!
       * @param grid points to the boundary's grid
       */
      explicit NorthBoundary(ConstGridPatchBasePointer grid) :
         Base(grid)
      {
      }

      //! Destruct a NorthBoundary
      virtual ~NorthBoundary()
      {
      }

      BOUNDARY_VERTEX_FUNCTIONS(ColumnIndex)
      BOUNDARY_EDGE_FUNCTIONS(ColumnIndex, (-1))
      BOUNDARY_VIRTUAL_FUNCTIONS

      //! Convert a ColumnIndex to a BoundaryIndex
      BoundaryIndex AsBoundaryIndex(ColumnIndex i) const
      {
         return this->ReverseIndex(BoundaryIndex(i.Get()));
      }

   protected:
      //! Converte a BoundaryIndex to a ColumnIndex
      ColumnIndex AsGridIndex(BoundaryIndex i) const
      {
         return ColumnIndex((this->ReverseIndex(i)).Get());
      }
   };

   //! Represent a GridRegion's EastBoundary
   /*!
    * @tparam Analysis is the analysis type
    */
   template<typename Analysis>
   struct EastBoundary: Boundary<Analysis> {
   private:
      BOUNDARY_USING

      public:
      //! Initialization
      /*!
       * @param grid points to the GridPatch this boundary belongs to.
       */
      explicit EastBoundary(ConstGridPatchBasePointer grid)
         : Base(grid)
      {
      }

      //! Destruct an EastBoundary
      virtual ~EastBoundary()
      {
      }

      BOUNDARY_VERTEX_FUNCTIONS(RowIndex)
      BOUNDARY_EDGE_FUNCTIONS(RowIndex, (0))
      BOUNDARY_VIRTUAL_FUNCTIONS

      //! Convert a ColumnIndex to a BoundaryIndex
      BoundaryIndex AsBoundaryIndex(RowIndex i) const {
         return BoundaryIndex(i.Get());
      }

   protected:
      //! Converte a BoundaryIndex to a ColumnIndex
      RowIndex AsGridIndex(BoundaryIndex i) const {
         return RowIndex(i.Get());
      }
   };

   //! Represent a GridRegion's SouthBoundary
   /*!
    * @tparam Analysis is the analysis type
    */
   template<typename Analysis>
   struct SouthBoundary: Boundary<Analysis> {
   private:
      BOUNDARY_USING

      public:
      //! Initialization
      /*!
       * @param grid points to the boundary's grid
       */
      explicit SouthBoundary(ConstGridPatchBasePointer grid)
         : Base(grid)
      {
      }

      //! Destruction
      virtual ~SouthBoundary()
      {
      }

      BOUNDARY_VERTEX_FUNCTIONS(ColumnIndex)
      BOUNDARY_EDGE_FUNCTIONS(ColumnIndex, (0))
      BOUNDARY_VIRTUAL_FUNCTIONS

      //! Convert a ColumnIndex to a BoundaryIndex
      BoundaryIndex AsBoundaryIndex(ColumnIndex i)  const {
         return BoundaryIndex(i.Get());
      }

   protected:
      //! Converte a BoundaryIndex to a ColumnIndex
      ColumnIndex AsGridIndex(BoundaryIndex i) const {
         return ColumnIndex(i.Get());
      }
   };

   //! Represent a GridRegion's WestBoundary
   /*!
    * @tparam Analysis is the analysis type
    */
   template<typename Analysis>
   struct WestBoundary: Boundary<Analysis> {
   private:
      BOUNDARY_USING

      public:
      //! Initialization
      /*!
       * @param grid points to the GridPatch this boundary belongs to.
       */
      explicit WestBoundary(ConstGridPatchBasePointer grid)
         : Base(grid)
      {
      }

      //! Destruction
      virtual ~WestBoundary()
      {
      }

      BOUNDARY_VERTEX_FUNCTIONS(RowIndex)
      BOUNDARY_EDGE_FUNCTIONS(RowIndex, (-1))
      BOUNDARY_VIRTUAL_FUNCTIONS

      //! Convert a ColumnIndex to a BoundaryIndex
      BoundaryIndex AsBoundaryIndex(RowIndex i) const {
         return this->ReverseIndex(BoundaryIndex(i.Get()));
      }

   protected:
      //! Converte a BoundaryIndex to a ColumnIndex
      RowIndex AsGridIndex(BoundaryIndex i) const {
         return RowIndex((this->ReverseIndex(i)).Get ());
      }
   };
}

#ifdef BOUNDARY_VERTEX_FUNCTIONS
#undef BOUNDARY_VERTEX_FUNCTIONS
#endif
#ifdef BOUNDARY_EDGE_FUNCTIONS
#undef BOUNDARY_EDGE_FUNCTIONS
#endif
#ifdef BOUNDARY_VIRTUAL_FUNCTIONS
#undef BOUNDARY_VIRTUAL_FUNCTIONS
#endif
#ifdef BOUNDARY_USING
#undef BOUNDARY_USING
#endif

//--------------------------------------------------------------
//
//   Function implementations
//
//--------------------------------------------------------------

#include "grid_patch_base.hpp"

namespace PACKAGE_NS {
   template<typename Analysis>
   inline void Boundary<Analysis>::CheckYourPartner() const
   {
      if (this->HasBoundaryCoupling())
         if (this->NumberOfVertices()
             != this->partner_boundary()->NumberOfVertices()) {
            boost::format error_message(
               "%1%: Attempt to couple edges having different number of vertices! First edge: %2% vertices; second edge: %3% vertices.");
#ifdef __GNUC__
#define FUNCTION_NAME __PRETTY_FUNCTION__
#else
#define FUNCTION_NAME (typeid(*this).name())
#endif
            throw RuntimeError(
               (error_message % FUNCTION_NAME % this->NumberOfVertices()
                % this->partner_boundary_->NumberOfVertices()).str());
#undef FUNCTION_NAME
         }
   }

   template<typename Analysis>
   inline BoundaryIndex Boundary<Analysis>::ReverseIndex(BoundaryIndex i) const
   {
      return BoundaryIndex(this->NumberOfVertices() - 1 - i.Get());
   }

   template<typename Analysis>
   inline void Boundary<Analysis>::SetBoundaryCoupling(
      const Boundary<Analysis>& other)
   {
      this->DeleteBoundaryCondition();
      partner_boundary_ = &other;
   }

   template<typename Analysis>
   inline void Boundary<Analysis>::DeleteBoundaryCoupling()
   {
      partner_boundary_ = ConstBoundaryPointer(nullptr);
   }

   template<typename Analysis>
   inline bool Boundary<Analysis>::HasBoundaryCoupling() const
   {
      return nullptr != partner_boundary_;
   }

   template<typename Analysis>
   inline void Boundary<Analysis>::SetBoundaryCondition(
      ConstBoundaryConditionPointer<Analysis> condition)
   {
      this->DeleteBoundaryCoupling();
      boundary_condition_ = condition;
   }

   template<typename Analysis>
   inline void Boundary<Analysis>::DeleteBoundaryCondition()
   {
      boundary_condition_ = ConstBoundaryConditionPointer<Analysis>(nullptr);
   }

   template<typename Analysis>
   inline bool Boundary<Analysis>::HasBoundaryCondition() const
   {
      return nullptr != boundary_condition_;
   }

   template <typename Analysis>
   inline auto Boundary<Analysis>::GetLoad(
      BoundaryIndex i) const -> LoadDefPointer
   {
      using std::get;

      if (this->HasBoundaryCoupling())
         return this->partner_boundary()->GetLoad(ReverseIndex(i));

      const auto row_and_column = AsGridVertexID(i);

      return this->grid_patch()->GetLoadHere(
         get<0>(row_and_column), get<1>(row_and_column));
   }

   template <typename Analysis>
   inline ConstMaterialPointer Boundary<Analysis>::GetMaterialAtVertex(
      BoundaryIndex i) const
   {
      using std::get;

      if (this->HasBoundaryCoupling())
         return this->partner_boundary()->GetMaterialAtVertex(
            ReverseIndex(i));

      const auto row_and_column = AsGridVertexID(i);

      return this->grid_patch()->GetMaterialHere(
         get<0>(row_and_column), get<1>(row_and_column));
   }

   template <typename Analysis>
   inline ConstMaterialPointer Boundary<Analysis>::GetMaterialAtEdge(
      BoundaryIndex i) const
   {
      using std::get;

      if (this->HasBoundaryCoupling())
         return this->partner_boundary()->GetMaterialAtEdge(
            ReverseIndex(i)-1);

      const auto edge_id = AsGridEdgeID(i);

      return this->grid_patch()->GetMaterialHere(
         get<0>(edge_id), get<1>(edge_id), get<2>(edge_id));
   }

   template<typename Analysis>
   inline VertexIndex Boundary<Analysis>::GlobalIndex(
      BoundaryIndex i, VertexIndexSetPointer previous_vertices) const
   {
      using std::get;

      if (this->HasBoundaryCoupling()) {
         const auto row_and_column_partner
            = this->partner_boundary()->AsGridVertexID(
               this->ReverseIndex(i));

         return this->partner_boundary()->grid_patch()->GlobalIndex(
            get<0>(row_and_column_partner),
            get<1>(row_and_column_partner),
            previous_vertices);
      }

      const auto row_and_column = this->AsGridVertexID(i);

      return this->grid_patch()->GlobalIndexHere(
         get<0>(row_and_column), get<1>(row_and_column));
   }

   template <typename Analysis>
   inline Length Boundary<Analysis>::AverageEdgeLengthPerpendicular(
      BoundaryIndex i) const
   {
      return this->HasBoundaryCoupling() ?
         (this->EdgeLengthPerpendicular(i)
          + partner_boundary()->EdgeLengthPerpendicular(ReverseIndex(i)))/2. :
         this->EdgeLengthPerpendicular(i);
   }

   template <typename Analysis>
   inline Area Boundary<Analysis>::CellSize(BoundaryIndex i) const
   {
      using std::get;

      if (this->HasBoundaryCoupling()) {
         return this->CellSizeInnerHalf(i)
            + this->partner_boundary()->CellSizeInnerHalf(
               this->ReverseIndex(i));
      }

      const auto grid_vertex_id = AsGridVertexID(i);

      return this->grid_patch()->CellSizeHere(
         get<0>(grid_vertex_id), get<1>(grid_vertex_id));
   }

   template <typename Analysis>
   inline Length Boundary<Analysis>::CrossingEdgeLength(BoundaryIndex i) const
   {
      using std::get;

      const auto edge_id = this->AsGridEdgeID(i);
      const Length crossing_edge_len_here
         = this->grid_patch()->CrossingEdgeLengthHere(
            get<0>(edge_id), get<1>(edge_id), get<2>(edge_id));

      if (this->HasBoundaryCoupling()) {
         const Length crossing_edge_len_partner
            = this->partner_boundary()->CrossingEdgeLength(
               this->ReverseIndex(i-1));

         return (crossing_edge_len_here + crossing_edge_len_partner)/2.;
      }

      return crossing_edge_len_here;
   }

   template <typename Analysis>
   CartesianSystem::Point Boundary<Analysis>::VertexPositionGlobal(
      BoundaryIndex i) const
   {
      using std::get;

      auto row_and_column = this->AsGridVertexID(i);
      return this->grid_patch()->VertexPositionGlobal(
         get<0>(row_and_column), get<1>(row_and_column));
   }

   //--------------------------------------------------------------
   //   North boundary
   //--------------------------------------------------------------

   template<typename Analysis>
   inline VertexIndex NorthBoundary<Analysis>::NumberOfVertices() const
   {
      return this->grid_patch()->NumberOfVertexColumns();
   }

   template<typename Analysis>
   inline VertexIndex NorthBoundary<Analysis>::NumberOfEdges() const
   {
      return this->grid_patch()->NumberOfVertexColumns() - 1;
   }

   template<typename Analysis>
   inline Length NorthBoundary<Analysis>::EdgeLengthPerpendicular(
      BoundaryIndex i) const
   {
      return this->grid_patch()->EdgeLength(
         this->grid_patch()->MaxRowIndex() - 1,
         this->AsGridIndex(i), Axis::kSecond);
   }

   template <typename Analysis>
   inline auto NorthBoundary<Analysis>::AsGridVertexID(
      BoundaryIndex i) const -> VertexID
   {
      return std::make_tuple(
         this->grid_patch()->MaxRowIndex(),
         this->AsGridIndex(i));
   }

   template <typename Analysis>
   inline auto NorthBoundary<Analysis>::AsGridEdgeID(
      BoundaryIndex i) const -> EdgeID
   {
      return std::make_tuple(
         this->grid_patch()->MaxRowIndex(),
         this->AsGridIndex(i-1), Axis::kFirst);
   }

   template <typename Analysis>
   inline Area NorthBoundary<Analysis>::CellSizeInnerHalf(
      BoundaryIndex i) const
   {
      using std::get;

      const auto grid_vertex_id = this->AsGridVertexID(i);

      return this->grid_patch()->CellSizeBottomHalf(
         get<0>(grid_vertex_id), get<1>(grid_vertex_id));
   }

   template <typename Analysis>
   Parallelity NorthBoundary<Analysis>::parallelity() const {
      return Parallelity::kParallel;
   }

   //--------------------------------------------------------------
   //   East boundary
   //--------------------------------------------------------------

   template<typename Analysis>
   inline VertexIndex EastBoundary<Analysis>::NumberOfVertices() const
   {
      return this->grid_patch()->NumberOfVertexRows();
   }

   template<typename Analysis>
   inline VertexIndex EastBoundary<Analysis>::NumberOfEdges() const
   {
      return this->grid_patch()->NumberOfVertexRows() - 1;
   }

   template<typename Analysis>
   inline Length EastBoundary<Analysis>::EdgeLengthPerpendicular(
      BoundaryIndex i) const
   {
      return this->grid_patch()->EdgeLength(
         AsGridIndex(i), this->grid_patch()->MaxColumnIndex() - 1,
         Axis::kFirst);
   }

   template <typename Analysis>
   inline auto EastBoundary<Analysis>::AsGridVertexID(
      BoundaryIndex i) const -> VertexID
   {
      return std::make_tuple(
         AsGridIndex(i),
         this->grid_patch()->MaxColumnIndex());
   }

   template <typename Analysis>
   inline auto EastBoundary<Analysis>::AsGridEdgeID(
      BoundaryIndex i) const -> EdgeID
   {
      return std::make_tuple(
         AsGridIndex(i),
         this->grid_patch()->MaxColumnIndex(),
         Axis::kSecond);
   }

   template <typename Analysis>
   inline Area EastBoundary<Analysis>::CellSizeInnerHalf(
      BoundaryIndex i) const
   {
      using std::get;

      const auto grid_vertex_id = this->AsGridVertexID(i);

      return this->grid_patch()->CellSizeLeftHalf(
         get<0>(grid_vertex_id), get<1>(grid_vertex_id));
   }

   template <typename Analysis>
   Parallelity EastBoundary<Analysis>::parallelity() const {
      return Parallelity::kParallel;
   }

   //--------------------------------------------------------------
   //   South boundary
   //--------------------------------------------------------------

   template<typename Analysis>
   inline VertexIndex SouthBoundary<Analysis>::NumberOfVertices() const
   {
      return this->grid_patch()->NumberOfVertexColumns();
   }

   template<typename Analysis>
   inline VertexIndex SouthBoundary<Analysis>::NumberOfEdges() const
   {
      return this->grid_patch()->NumberOfVertexColumns() - 1;
   }

   template<typename Analysis>
   inline Length SouthBoundary<Analysis>::EdgeLengthPerpendicular(
      BoundaryIndex i) const
   {
      return this->grid_patch()->EdgeLength(
         RowIndex(0), this->AsGridIndex(i), Axis::kSecond);
   }

   template <typename Analysis>
   inline auto SouthBoundary<Analysis>::AsGridVertexID(
      BoundaryIndex i) const -> VertexID
   {
      return std::make_tuple(RowIndex(0), AsGridIndex(i));
   }

   template <typename Analysis>
   inline auto SouthBoundary<Analysis>::AsGridEdgeID(
      BoundaryIndex i) const -> EdgeID
   {
      return std::make_tuple(RowIndex(0), AsGridIndex(i), Axis::kFirst);
   }

   template <typename Analysis>
   inline Area SouthBoundary<Analysis>::CellSizeInnerHalf(
      BoundaryIndex i) const
   {
      using std::get;

      const auto grid_vertex_id = this->AsGridVertexID(i);

      return this->grid_patch()->CellSizeTopHalf(
         get<0>(grid_vertex_id), get<1>(grid_vertex_id));
   }

   template <typename Analysis>
   Parallelity SouthBoundary<Analysis>::parallelity() const {
      return Parallelity::kAntiparallel;
   }

   //--------------------------------------------------------------
   //   West boundary
   //--------------------------------------------------------------

   template<typename Analysis>
   inline VertexIndex WestBoundary<Analysis>::NumberOfVertices() const
   {
      return this->grid_patch()->NumberOfVertexRows();
   }

   template<typename Analysis>
   inline VertexIndex WestBoundary<Analysis>::NumberOfEdges() const
   {
      return this->grid_patch()->NumberOfVertexRows() - 1;
   }

   template<typename Analysis>
   inline Length WestBoundary<Analysis>::EdgeLengthPerpendicular(
      BoundaryIndex i) const
   {
      return this->grid_patch()->EdgeLength(
         AsGridIndex(i), ColumnIndex(0), Axis::kFirst);
   }

   template <typename Analysis>
   inline auto WestBoundary<Analysis>::AsGridVertexID(
      BoundaryIndex i) const -> VertexID
   {
      return std::make_tuple(this->AsGridIndex(i), ColumnIndex(0));
   }

   template <typename Analysis>
   inline auto WestBoundary<Analysis>::AsGridEdgeID(
      BoundaryIndex i) const -> EdgeID
   {
      return std::make_tuple(
         this->AsGridIndex(i-1), ColumnIndex(0), Axis::kSecond);
   }

   template <typename Analysis>
   inline Area WestBoundary<Analysis>::CellSizeInnerHalf(
      BoundaryIndex i) const
   {
      using std::get;

      const auto grid_vertex_id = this->AsGridVertexID(i);

      return this->grid_patch()->CellSizeRightHalf(
         get<0>(grid_vertex_id), get<1>(grid_vertex_id));
   }

   template <typename Analysis>
   Parallelity WestBoundary<Analysis>::parallelity() const {
      return Parallelity::kAntiparallel;
   }

#ifndef INSTANTIATE_
#define EXTERN_TEMPLATE_EDGES(NEUMANN_ANALYSIS)             \
   extern template struct EastBoundary<NEUMANN_ANALYSIS>;   \
   extern template struct NorthBoundary<NEUMANN_ANALYSIS>;  \
   extern template struct WestBoundary<NEUMANN_ANALYSIS>;   \
   extern template struct SouthBoundary<NEUMANN_ANALYSIS>;

   EXTERN_TEMPLATE_EDGES(TimeTransientAnalysis<>)
   EXTERN_TEMPLATE_EDGES(MagnetoHarmonicAnalysis<>)
#undef EXTERN_TEMPLATE_EDGES
#endif

}

#endif // BOUNDARY_HPP_FGr2WLZN_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
