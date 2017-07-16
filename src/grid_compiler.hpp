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
 *  \brief grid compiler
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef GRID_COMPILER_HPP_pZgso9bJ_
#define GRID_COMPILER_HPP_pZgso9bJ_

#include "config.h"

#include "boundary_vertex_factory.hpp"
#include "grids-fwd.hpp"
#include "grid_graph-fwd.hpp"
#include "loki_visitor_extra.hpp"

#define GRID_ENTITY_CREATOR_TYPES                                       \
   using BoundaryVertexFactoryType = BoundaryVertexFactory<Analysis>;	\
   using EdgePropertyIndexVector = std::vector<std::size_t>;            \
   using EdgePropertyVector = std::vector<EdgeProps>;                   \
   using GridEdgeDescriptor = std::pair<VertexIndex, VertexIndex>;		\
   using EdgeVector = std::vector<GridEdgeDescriptor>;                  \
   using GridGraphPtr = GridGraphPointer<Analysis>;                     \
   using GridGraphType = GridGraph<Analysis>;                           \
   using StandardVertexType = StandardVertex<Analysis>;                 \
   using VertexPtr = VertexPointer<Analysis>;                           \
   using VertexPropertyVector = std::vector<VertexPtr>;                 \
   using size_t = std::size_t;

#define GRID_ENTITY_CREATOR_FUNCTION_DECLS                              \
   template <typename CoordinateSystem>                                 \
   void CollectEdge(                                                    \
      const GridPatch<CoordinateSystem, Analysis>& grid_patch,          \
      RowIndex row_index, ColumnIndex column_index, Axis axis,          \
      EdgeVector& edges, EdgePropertyIndexVector& edge_property_index,	\
      EdgePropertyVector& edge_properties) const;                       \
   template <typename CoordinateSystem>                                 \
   void CreateInnerVertex(                                              \
      const GridPatch<CoordinateSystem, Analysis>& grid_patch,          \
      RowIndex row_index, ColumnIndex column_index,                     \
      VertexPtr& new_vertex) const;

namespace PACKAGE_NS {
   namespace fi {
      namespace detail {
         //! Methods for creating edge and vertex objects
         /*!
          * \tparam Analysis is the analysis type
          */
         template <typename Analysis>
         struct GridEntityCreation
         {
         protected:
            using EdgeProps = GridEdgePointer<Analysis>;
            GRID_ENTITY_CREATOR_TYPES

            public:
            GRID_ENTITY_CREATOR_FUNCTION_DECLS

            //! Calculate the required size for the EdgePropertyVector
            /*!
             * @param number_of_edges number of edges
             */
            size_t EdgePropertyVectorRequiredSize(
               size_t number_of_edges) {
               return number_of_edges;
            }
         };
      }
   }

   namespace detail {
      //! Compile a grid graph from a grid setup
      /*!
       * \tparam Analysis is the analysis type
       *
       * \tparam PlugIn defines how grid entities (vertices and edges)
       * are created.
       */
      template <typename Analysis, template <typename> class PlugIn>
      struct GridCompiler :
         PlugIn<Analysis>,
         Loki::BaseVisitor,
         Loki::Visitor<GridAssembly, void, true>,
         Loki::Visitor<GridPatch<CartesianSystem, Analysis>, void, true>,
         Loki::Visitor<GridPatch<PolarSystem, Analysis>, void, true>
      {
      private:
         using EntityCreator = PlugIn<Analysis>;
         using BoundaryVertexFactory =
            typename EntityCreator::BoundaryVertexFactoryType;
         //! Type for describing a grid edge
         //using GridEdgeDescriptor = typename EntityCreator::GridEdgeDescriptor;
         //! Vector for storing grid edges
         using EdgeVector = typename EntityCreator::EdgeVector;
         using EdgeProps = typename EntityCreator::EdgeProps;
         using VertexPtr = typename EntityCreator::VertexPtr;
         //! Vector for storing edge properties
         using EdgePropertyVector =
            typename EntityCreator::EdgePropertyVector;
         //! Vector for storing vertex properties
         using VertexPropertyVector =
            typename EntityCreator::VertexPropertyVector;
         using BoundaryType = Boundary<Analysis>;
         using EdgePropertyIndexVector =
            typename EntityCreator::EdgePropertyIndexVector;
         using Time = typename Analysis::Time;
         using GridPtr = ConstGridPointer;
         using GridGraphPtr = typename EntityCreator::GridGraphPtr;
         using GridGraphType = typename EntityCreator::GridGraphType;
         using StandardVertexType =
            typename EntityCreator::StandardVertexType;

      public:
         //! Initialization
         /*!
          * @param grid points to some grid (e.g. GridPatch or
          * GridAssembly)
          */
         explicit GridCompiler(GridPtr grid)
            : grid_(grid), edges_(nullptr), edge_property_index_(nullptr),
              edge_properties_(nullptr), vertex_properties_(nullptr)
         {
         }

         GridCompiler(const GridCompiler& other) = delete;

         //! Visit a GridAssembly
         void Visit(const GridAssembly&) override;

         //! Visit a GridPatch
         void Visit(
            const GridPatch<CartesianSystem, Analysis>& grid_patch) override
         {
            DoVisit(grid_patch);
         }

         //! Visit a GridPatch
         void Visit(
            const GridPatch<PolarSystem, Analysis>& grid_patch) override
         {
            DoVisit(grid_patch);
         }

         //! Visit a GridPatch
         template <typename CoordinateSystem>
         void DoVisit(const GridPatch<CoordinateSystem, Analysis>&);

         //! Build a GridGraph
         GridGraphPtr operator()(Time);

         //! Build a GridGraph
         GridGraphPtr operator()();

      private:
         template<typename CoordinateSystem>
         void CollectVertices(
            const GridPatch<CoordinateSystem, Analysis>& grid_patch);

         template <typename CoordinateSystem>
         inline void InsertBoundaryVertex(
            const GridPatch<CoordinateSystem, Analysis>& grid_patch,
            VertexPtr& vertex, RowIndex row_index, ColumnIndex column_index);

         //! Create a vertex object for a GridPatch'es corner vertex
         /*!
          * Determine the correct vertex type for a GridPatch'es corner vertex:
          *
          * 1. If the vertex is coupled to another vertex, it can be ignored
          * and nullptr is returned.
          *
          * 2. If one of the boundaries has a boundary condition set,
          * this condition is used to create the vertex.
          *
          * 3. If both boundaries have a condition set, the conditions
          * must "fight" and the stronger condition is used to create the
          * vertex.
          *
          * 4. If no condition is set, a StandardVertex is used.
          *
          * @param in_boundary is the boundary having the maximum vertex
          * index at the corner
          *
          * @param out_boundary is the boundary with vertex index zero at
          * the corner
          *
          * @param new_vertex[out] An appropriate vertex object or
          * nullptr if the corner vertex is coupled to another vertex
          */
         inline void CreateVertex4GridPatchCorner(
            const BoundaryType& in_boundary,
            const BoundaryType& out_boundary,
            BoundaryVertexFactory& new_vertex) const;

         template<typename CoordinateSystem>
         void CollectEdges(
            const GridPatch<CoordinateSystem, Analysis>& grid_patch);

         //! Remove "holes" in the vertex numbering
         /*!
          * @tparam Analysis is the analysis type
          */
         inline void CreateDenseVertexNumbering();

         //! Figure out a good equation numbering scheme
         /*!
          * Figures out a good equation numbering scheme and sets the
          * equation number in the vertex object.
          *
          * @tparam Analysis is the Analysis type
          * @param[in] edges is a vector with all edges
          * @param vertices is a vector with all vertices
          * @return the number of equations
          */
         inline EquationNumber FindEquationNumbering();

      private:
         GridPtr grid_;
         EdgeVector* edges_;
         EdgePropertyIndexVector* edge_property_index_;
         EdgePropertyVector* edge_properties_;
         VertexPropertyVector* vertex_properties_;
      };
   }
}

//--------------------------------------------------------------
//
//   Function implementations
//
//--------------------------------------------------------------

#include <algorithm>
#include <map>
#include <vector>
#include <utility>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>

#include "edge_maker.hpp"
#include "grid.hpp"
#include "grid_assembly.hpp"
#include "grid_graph.hpp"
#include "grid_patch.hpp"
#include "standard_vertex.hpp"

namespace PACKAGE_NS {
   namespace fi {
      namespace detail {
         template <typename Analysis>
         template <typename CoordinateSystem>
         inline void GridEntityCreation<Analysis>::CollectEdge(
            const GridPatch<CoordinateSystem, Analysis>& grid_patch,
            RowIndex row_index, ColumnIndex column_index, Axis axis,
            EdgeVector& edges, EdgePropertyIndexVector& edge_property_index,
            EdgePropertyVector& edge_properties) const
         {
            const VertexIndex start_vertex
               = grid_patch.GlobalIndex(row_index, column_index);
            const VertexIndex end_vertex
               = grid_patch.GlobalIndex(
                  row_index + (Axis::kSecond == axis),
                  column_index + (Axis::kFirst == axis));

            edges.emplace_back(start_vertex, end_vertex);
            edges.emplace_back(end_vertex, start_vertex);

            EdgeMaker<Analysis> edge_maker(row_index, column_index, axis);

            edge_properties.push_back(edge_maker.Visit(grid_patch));

            std::size_t ep_index = edge_properties.size() - 1;

            edge_property_index.emplace_back(ep_index);
            edge_property_index.emplace_back(ep_index);
         }

         template <typename Analysis>
         template <typename CoordinateSystem>
         inline void
         GridEntityCreation<Analysis>::CreateInnerVertex(
            const GridPatch<CoordinateSystem, Analysis>& grid_patch,
            RowIndex row_index, ColumnIndex column_index,
            VertexPtr& new_vertex) const
         {
            new_vertex = std::make_unique<StandardVertex<Analysis> >(
               grid_patch.VertexPositionGlobal(row_index, column_index),
               grid_patch.GetLoad(row_index, column_index),
               grid_patch.CellSize(row_index, column_index));
         }
      }
   }

   namespace detail {
      template <typename Analysis, template <typename> class PlugIn>
      void GridCompiler<Analysis, PlugIn>::Visit(
         const GridAssembly& grid_assembly)
      {
         for (auto iter = grid_assembly.grids().begin(),
                 end = grid_assembly.grids().end();
              iter != end; ++iter) {
            (*iter)->Accept(*this);
         }
      }

      template <typename Analysis, template <typename> class PlugIn>
      template <typename CoordinateSystem>
      void GridCompiler<Analysis, PlugIn>::DoVisit(
         const GridPatch<CoordinateSystem, Analysis>& grid_patch)
      {
         grid_patch.CheckBoundaryCompability();
         this->CollectEdges(grid_patch);
         this->CollectVertices(grid_patch);
      }

      template <typename Analysis, template <typename> class PlugIn>
      inline auto GridCompiler<Analysis, PlugIn>::operator()() -> GridGraphPtr
      {
         using Time = typename Analysis::Time;

         return this->operator()(Time::from_value(1.));
      }

      template <typename Analysis, template <typename> class PlugIn>
      auto GridCompiler<Analysis, PlugIn>::operator()(
         typename Analysis::Time time_step) -> GridGraphPtr
      {
         VertexPropertyVector vertex_properties(
            grid_->NumberOfVerticesBrutto());
         EdgeVector edges;
         EdgePropertyVector edge_properties;
         EdgePropertyIndexVector edge_property_index;
         this->vertex_properties_ = &vertex_properties;
         this->edges_ = &edges;
         this->edge_properties_ = &edge_properties;
         this->edge_property_index_ = &edge_property_index;
         const auto number_of_edges_netto = grid_->NumberOfEdgesNetto();
         const auto number_of_edge_elements =
            this->EdgePropertyVectorRequiredSize(number_of_edges_netto);
         edges_->reserve(2*number_of_edges_netto);
         edge_property_index_->reserve(2*number_of_edges_netto);
         edge_properties_->reserve(number_of_edge_elements);
         grid_->SetFirstVertexIndex();
         grid_->Accept(*this);
         this->CreateDenseVertexNumbering();
         const auto number_of_equations = this->FindEquationNumbering();

         // Step 4: Create the grid graph
         auto grid_graph = std::make_shared<GridGraphType>(
            edges_->begin(), edges_->end(), edge_property_index_->begin(),
            edge_properties_, vertex_properties_->size(),
            GraphProperties<Analysis>(number_of_equations, time_step));

         // Attach vertex properties to the grid graph
         for (auto v = vertices(*grid_graph); v.first != v.second; ++(v.first))
         {
#ifdef NDEBUG
            std::swap((*grid_graph)[*(v.first)],
                      (*vertex_properties_)[*(v.first)]);
#else
            std::swap((*grid_graph)[*(v.first)],
                      vertex_properties_->at(*(v.first)));
#endif
         }

         return grid_graph;
      }

      template <typename Analysis, template <typename> class PlugIn>
      template <typename CoordinateSystem>
      void GridCompiler<Analysis, PlugIn>::CollectVertices(
         const GridPatch<CoordinateSystem, Analysis>& grid_patch)
      {
         BoundaryVertexFactory boundary_vertex_factory;

         const auto max_row_index = grid_patch.MaxRowIndex();
         const auto max_column_index = grid_patch.MaxColumnIndex();

         // *** Corners ***
         //        South-West
         CreateVertex4GridPatchCorner(
            grid_patch.west_boundary(), grid_patch.south_boundary(),
            boundary_vertex_factory);
         InsertBoundaryVertex(
            grid_patch, boundary_vertex_factory.new_vertex_object(),
            RowIndex(0), ColumnIndex(0));
         //        South-East
         CreateVertex4GridPatchCorner(
            grid_patch.south_boundary(), grid_patch.east_boundary(),
            boundary_vertex_factory);
         InsertBoundaryVertex(
            grid_patch, boundary_vertex_factory.new_vertex_object(),
            RowIndex(0), max_column_index);
         //        North-East
         CreateVertex4GridPatchCorner(
            grid_patch.east_boundary(), grid_patch.north_boundary(),
            boundary_vertex_factory);
         InsertBoundaryVertex(
            grid_patch, boundary_vertex_factory.new_vertex_object(),
            max_row_index, max_column_index);
         //        North-West
         CreateVertex4GridPatchCorner(
            grid_patch.north_boundary(), grid_patch.west_boundary(),
            boundary_vertex_factory);
         InsertBoundaryVertex(
            grid_patch, boundary_vertex_factory.new_vertex_object(),
            max_row_index, ColumnIndex(0));

         // *** East and west boundary ***
         for (RowIndex row_index(1); row_index < max_row_index; ++row_index)
         {
            // east boundary
            boundary_vertex_factory(
               grid_patch.east_boundary(), row_index);
            InsertBoundaryVertex(
               grid_patch, boundary_vertex_factory.new_vertex_object(),
               row_index, max_column_index);

            // west boundary
            boundary_vertex_factory(grid_patch.west_boundary(), row_index);
            InsertBoundaryVertex(
               grid_patch, boundary_vertex_factory.new_vertex_object(),
               row_index, ColumnIndex(0));
         }

         // *** North and south boundary ***
         for (ColumnIndex column_index(1);
              column_index < max_column_index; ++column_index)
         {
            // north boundary
            boundary_vertex_factory(
               grid_patch.north_boundary(), column_index);
            InsertBoundaryVertex(
               grid_patch, boundary_vertex_factory.new_vertex_object(),
               max_row_index, column_index);
            // south boundary
            boundary_vertex_factory(
               grid_patch.south_boundary(), column_index);
            InsertBoundaryVertex(
               grid_patch, boundary_vertex_factory.new_vertex_object(),
               RowIndex(0), column_index);
         }

         // Inner vertex rows
         for (RowIndex row_index(1);
              row_index < max_row_index; ++row_index)
         {
            // inner vertices
            for (ColumnIndex column_index(1);
                 column_index < max_column_index; ++column_index)
            {
               VertexPtr new_vertex(nullptr);
               this->CreateInnerVertex(
                  grid_patch, row_index, column_index, new_vertex);
               auto global_vertex_index = grid_patch.GlobalIndex(
                  row_index, column_index);
#ifdef NDEBUG
               std::swap((*vertex_properties_)[global_vertex_index],
                         new_vertex);
#else
               std::swap((*vertex_properties_).at(global_vertex_index),
                         new_vertex);
#endif
            }
         }
      }

      template <typename Analysis, template <typename> class PlugIn>
      template <typename CoordinateSystem>
      inline void GridCompiler<Analysis, PlugIn>::InsertBoundaryVertex(
         const GridPatch<CoordinateSystem, Analysis>& grid_patch,
         VertexPtr& vertex, RowIndex row_index, ColumnIndex column_index)
      {
         if (nullptr != vertex)
         {
            const auto global_vertex_index = grid_patch.GlobalIndex(
               row_index, column_index);
#ifdef NDEBUG
            std::swap((*vertex_properties_)[global_vertex_index], vertex);
#else
            std::swap((*vertex_properties_).at(global_vertex_index), vertex);
#endif
         }
      }

      template<typename Analysis, template <typename> class PlugIn>
      inline void
      GridCompiler<Analysis, PlugIn>::CreateVertex4GridPatchCorner(
         const Boundary<Analysis>& in_boundary,
         const Boundary<Analysis>& out_boundary,
         BoundaryVertexFactory& boundary_vertex_factory) const
      {
         const bool in_has_boundary_coupling =
            in_boundary.HasBoundaryCoupling();
         const bool out_has_boundary_coupling =
            out_boundary.HasBoundaryCoupling();

         const bool in_has_boundary_condition =
            in_boundary.HasBoundaryCondition();
         const bool out_has_boundary_condition =
            out_boundary.HasBoundaryCondition();

         if (in_has_boundary_coupling || out_has_boundary_coupling) {
            VertexPtr new_vertex = nullptr;
            boundary_vertex_factory.set_new_vertex_object(new_vertex);
         }
         else if (in_has_boundary_condition and out_has_boundary_condition)
         {
            ConstBoundaryConditionPointer<Analysis>
               in_boundary_condition = in_boundary.boundary_condition(),
               out_boundary_condition = out_boundary.boundary_condition();

            const short
               in_priority = in_boundary_condition->priority(),
               out_priority = out_boundary_condition->priority();

            if (in_priority > out_priority)
               boundary_vertex_factory(
                  in_boundary, in_boundary.ReverseIndex(BoundaryIndex(0)));
            else
               boundary_vertex_factory(out_boundary, BoundaryIndex(0));
         }
         else if (in_has_boundary_condition)
            boundary_vertex_factory(
               in_boundary, in_boundary.ReverseIndex(BoundaryIndex(0)));
         else
            boundary_vertex_factory(out_boundary, BoundaryIndex(0));
      }

      template<typename Analysis, template <typename> class PlugIn>
      template<typename CoordinateSystem>
      void GridCompiler<Analysis, PlugIn>::CollectEdges(
         const GridPatch<CoordinateSystem, Analysis>& grid_patch)
      {
         const bool east_has_boundary_coupling
            = grid_patch.east_boundary().HasBoundaryCoupling();
         const bool north_has_boundary_coupling
            = grid_patch.north_boundary().HasBoundaryCoupling();
         const bool west_has_boundary_coupling
            = grid_patch.west_boundary().HasBoundaryCoupling();
         const bool south_has_boundary_coupling
            = grid_patch.south_boundary().HasBoundaryCoupling();

         const auto max_row_index = grid_patch.MaxRowIndex();
         const auto max_column_index = grid_patch.MaxColumnIndex();

         // South-west-corner
         if (not south_has_boundary_coupling)
            this->CollectEdge(
               grid_patch, RowIndex(0), ColumnIndex(0), Axis::kFirst,
               (*edges_), (*edge_property_index_), (*edge_properties_));
         if (not west_has_boundary_coupling)
            this->CollectEdge(
               grid_patch, RowIndex(0), ColumnIndex(0), Axis::kSecond,
               (*edges_), (*edge_property_index_), (*edge_properties_));

         // South edge
         for (ColumnIndex i(1); i < max_column_index; ++i) {
            if (not south_has_boundary_coupling)
               this->CollectEdge(
                  grid_patch, RowIndex(0), i, Axis::kFirst,
                  (*edges_), (*edge_property_index_), (*edge_properties_));
            this->CollectEdge(
               grid_patch, RowIndex(0), i, Axis::kSecond,
               (*edges_), (*edge_property_index_), (*edge_properties_));
         }

         // South-east corner
         if (not east_has_boundary_coupling)
            this->CollectEdge(
               grid_patch, RowIndex(0), max_column_index, Axis::kSecond,
               (*edges_), (*edge_property_index_), (*edge_properties_));

         // Inner vertex rows
         for (RowIndex row_index(1);
              row_index < max_row_index; ++row_index) {
            // West edge
            if (not west_has_boundary_coupling)
               this->CollectEdge(
                  grid_patch, row_index, ColumnIndex(0), Axis::kSecond,
                  (*edges_), (*edge_property_index_), (*edge_properties_));
            this->CollectEdge(
               grid_patch, row_index, ColumnIndex(0), Axis::kFirst,
               (*edges_), (*edge_property_index_), (*edge_properties_));

            for (ColumnIndex column_index(1);
                 column_index < max_column_index; ++column_index) {
               this->CollectEdge(
                  grid_patch, row_index, column_index, Axis::kFirst,
                  (*edges_), (*edge_property_index_), (*edge_properties_));
               this->CollectEdge(
                  grid_patch, row_index, column_index, Axis::kSecond,
                  (*edges_), (*edge_property_index_), (*edge_properties_));
            }

            // East edge
            if (not east_has_boundary_coupling)
               this->CollectEdge(
                  grid_patch, row_index, max_column_index, Axis::kSecond,
                  (*edges_), (*edge_property_index_), (*edge_properties_));
         }

         // north edge
         if (not north_has_boundary_coupling)
            for (ColumnIndex column_index(0);
                 column_index < max_column_index; ++column_index)
               this->CollectEdge(
                  grid_patch, max_row_index, column_index, Axis::kFirst,
                  (*edges_), (*edge_property_index_), (*edge_properties_));
      }

      template <typename Analysis, template <typename> class PlugIn>
      inline void GridCompiler<Analysis, PlugIn>::CreateDenseVertexNumbering()
      {
         std::map<VertexIndex, VertexIndex> sparse_to_dense_mapping;

         // Create mapping
         //
         // Because the graph is bidirectional it is sufficient to
         // iterate over the edge'es source vertices
         for (auto iter = edges_->begin(), end = edges_->end(); iter != end;
              ++iter) {
            if (not sparse_to_dense_mapping.count(iter->first)) {
               sparse_to_dense_mapping.emplace(iter->first, 0);
            }
         }

         VertexIndex dense_vertex_number_counter = 0;
         for (auto iter = sparse_to_dense_mapping.begin(), end =
                 sparse_to_dense_mapping.end(); iter != end; ++iter) {
            iter->second = dense_vertex_number_counter++;
         }

         // Renumber edges
         EdgeVector edges_new_numbering;
         edges_new_numbering.reserve(edges_->size());
         for (auto iter = edges_->begin(), end = edges_->end(); iter != end;
              ++iter) {
            edges_new_numbering.emplace_back(
               sparse_to_dense_mapping[iter->first],
               sparse_to_dense_mapping[iter->second]);
         }

         // Renumber/shift vertices
         VertexPropertyVector vertices_new_numbering(
            dense_vertex_number_counter);
         for (auto iter = sparse_to_dense_mapping.begin(), end =
                 sparse_to_dense_mapping.end(); iter != end; ++iter) {
#ifdef NDEBUG
            std::swap(vertices_new_numbering[iter->second],
                      (*vertex_properties_)[iter->first]);
#else
            std::swap(vertices_new_numbering.at(iter->second),
                      vertex_properties_->at(iter->first));
#endif
         }

         std::swap((*edges_), edges_new_numbering);
         std::swap((*vertex_properties_), vertices_new_numbering);
      }

      template<typename Analysis, template <typename> class PlugIn>
      inline EquationNumber
      GridCompiler<Analysis, PlugIn>::FindEquationNumbering()
      {
         using AdjacencyGraph =
            boost::adjacency_list<
               boost::vecS, boost::vecS, boost::bidirectionalS>;
         using IntVector = std::vector<int>;

         EdgeVector active_edges;
         active_edges.reserve(edges_->size());

         std::copy_if(
            edges_->begin(), edges_->end(),
            std::back_inserter(active_edges),
            [this](typename EdgeVector::value_type e) {
#ifdef NDEBUG
               return (
                  (nullptr != dynamic_cast<StandardVertexType* >(
                     ((*vertex_properties_)[e.first]).get()))
                  and (nullptr != dynamic_cast<StandardVertexType* >(
                          ((*vertex_properties_)[e.second]).get())));
#else
               return (
                  (nullptr != dynamic_cast<StandardVertexType*>(
                     (vertex_properties_->at(e.first)).get()))
                  and (nullptr != dynamic_cast<StandardVertexType*>(
                          (vertex_properties_->at(e.second)).get())));
#endif
            });

         std::map<VertexIndex, VertexIndex> sparse_to_dense_mapping;

         for (auto iter = active_edges.begin(), end = active_edges.end();
              iter != end; ++iter) {
            if (not sparse_to_dense_mapping.count(iter->first))
               sparse_to_dense_mapping.emplace(iter->first, 0);
         }

         VertexIndex dense_vertex_number_counter = 0;
         for (auto iter = sparse_to_dense_mapping.begin(), end =
                 sparse_to_dense_mapping.end(); iter != end; ++iter) {
            iter->second = dense_vertex_number_counter++;
         }

         AdjacencyGraph g;

         for (auto iter = active_edges.begin(), end = active_edges.end();
              iter != end; ++iter) {
            add_edge(sparse_to_dense_mapping[iter->first],
                     sparse_to_dense_mapping[iter->second], g);
         }

         const std::size_t number_of_vertices = num_vertices(g);

         IntVector inverse_perm(number_of_vertices, 0);
         IntVector perm(number_of_vertices, 0);
         IntVector supervertex_sizes(number_of_vertices, 1); // init has to be 1
         auto id = get(boost::vertex_index, g);
         IntVector outdegree(number_of_vertices, 0);
         constexpr int delta = 0;

         minimum_degree_ordering(
            g,
            make_iterator_property_map(outdegree.begin(), id),
            make_iterator_property_map(inverse_perm.begin(), id),
            make_iterator_property_map(perm.begin(), id),
            make_iterator_property_map(supervertex_sizes.begin(), id), delta,
            id);

         for (auto iter = sparse_to_dense_mapping.begin(), end =
                 sparse_to_dense_mapping.end(); iter != end; ++iter) {
            static_cast<StandardVertexType*>(
               (*vertex_properties_)[iter->first].get())
               ->set_equation_number(perm[iter->second]);
         }

         return number_of_vertices; // = number of equations
      }
   }

   namespace fi {
      template <typename Analysis>
      using GridCompiler = ::PACKAGE_NS::detail::GridCompiler<
         Analysis, detail::GridEntityCreation>;
   }

   namespace detail {
#ifndef INSTANTIATE_
      extern template struct GridCompiler<
         MagnetoHarmonicAnalysis<>,
         ::PACKAGE_NS::fi::detail::GridEntityCreation >;
      extern template struct GridCompiler<
         TimeTransientAnalysis<>,
         ::PACKAGE_NS::fi::detail::GridEntityCreation >;
#endif
   }
}

#ifdef GRID_ENTITY_CREATOR_TYPES
#undef GRID_ENTITY_CREATOR_TYPES
#endif
#ifdef GRID_ENTITY_CREATOR_FUNCTION_DECLS
#undef GRID_ENTITY_CREATOR_FUNCTION_DECLS
#endif

#endif // GRID_COMPILER_HPP_pZgso9bJ_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
