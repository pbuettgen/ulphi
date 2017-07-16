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
 *  \brief GridGraph and related stuff.
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef GRID_GRAPH_HPP_YV9h1g37_
#define GRID_GRAPH_HPP_YV9h1g37_

#include <functional>
#include <memory>

#include "config.h"

#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/graphviz.hpp>

#include "analysis.hpp"
#include "boost_iterator_extra.hpp"
#include "coordinate_systems.hpp"
#include "exception.hpp"
#include "types.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   //! (Global) properties attached to the graph
   /*!
    * @tparam Analysis is the analysis type
    *
    * @relates GridGraph
    */
   template<typename Analysis>
   struct GraphProperties {
   protected:
      using Time = typename Analysis::Time;

   public:
      //! Initialization
      /*
       * \param number_of_equations number of equations
       *
       * \param time_step time step (frequency step for a
       * magnetoharmonic analysis).
       */
      GraphProperties(
         EquationNumber number_of_equations, Time time_step)
         : number_of_equations_(number_of_equations),
           time_step_(time_step),
           time_step_counter_(0)
      {
      }

      //! Total number of equations
      /*!
       * @return total number of equations
       */
      EquationNumber number_of_equations() const
      {
         return number_of_equations_;
      }

      //! Set total number of equations
      /*!
       * @param number_of_equations is the new number of equations
       */
      void set_number_of_equations(EquationNumber number_of_equations)
      {
         number_of_equations_ = number_of_equations;
      }

      //! Time step / base frequency
      /*!
       * @return time step (time transient analysis) or base frequency
       * (time harmonic analysis)
       */
      const Time& time_step() const
      {
         return time_step_;
      }

      //! Get the current time
      Time current_time() const {
         return static_cast<typename Time::value_type>(
            time_step_counter_) * time_step_;
      }

      //! Procede with the next time step
      void NextTimeStep() {
         ++time_step_counter_;
      }

      //! Set time step / base frequency
      /*!
       * @param time_step the new time step / base frequency
       */
      void set_time_step(const Time& time_step)
      {
         time_step_ = time_step;
      }

   private:
      EquationNumber number_of_equations_;
      Time time_step_;
      std::size_t time_step_counter_;
   };

   namespace detail {
      //! Grid graph
      /*!
       * @tparam Analysis is the analysis type
       * @tparam VertexPointer pointer to a Vertex
       * @tparam EdgeProperties edge property type
       */
      template <typename Analysis,
                template <typename> class VertexPointer,
                template <typename> class EdgeProperties>
      struct GridGraph : boost::compressed_sparse_row_graph<
         boost::directedS,
         VertexPointer<Analysis>,
         std::size_t,
         GraphProperties<Analysis> >
      {
      private:
         using GraphPropertyType = GraphProperties<Analysis>;

         using Base = boost::compressed_sparse_row_graph<
            boost::directedS, VertexPointer<Analysis>,
            std::size_t, GraphPropertyType >;

         using EdgeProps = EdgeProperties<Analysis>;
         using EdgePropertyVector = std::vector<EdgeProps >;

         using graph_bundled = typename Base::graph_bundled;
         using VertexProperty = VertexPointer<Analysis>;

      public:
         //! Analysis type
         using AnalysisT = Analysis;

         //! edge descriptor
         using edge_descriptor = typename Base::edge_descriptor;

         using vertex_descriptor = typename Base::vertex_descriptor;

      public:
         //! Initialization
         /*!
          * @param edge_begin marks the start of a sequence with edges
          *
          * @param edge_end marks the end of the edge sequence
          *
          * @param edge_property_index_iterator sequence with indices
          * of the edge'es properties
          *
          * @param[in,out] edge_properties the graph takes over the
          * content of this vector, so it is empty after the call to
          * this constructor.
          *
          * @param numverts number of vertices
          *
          * @param graph_properties global properties
          */
         template<typename MultiPassInputIterator,
                  typename EdgePropertyIndexIterator>
         GridGraph(
            MultiPassInputIterator edge_begin, MultiPassInputIterator edge_end,
            EdgePropertyIndexIterator edge_property_index_iterator,
            EdgePropertyVector* edge_properties, std::size_t numverts,
            const GraphPropertyType& graph_properties)
            : Base(boost::edges_are_unsorted_multi_pass,
                   edge_begin, edge_end, edge_property_index_iterator,
                   numverts, graph_properties)
         {
            std::swap(edge_properties_, *edge_properties);
         }

         VertexProperty& operator[](const vertex_descriptor& vertex_id) {
            return Base::operator[](vertex_id);
         }

         const VertexProperty& operator[](
            const vertex_descriptor& vertex_id) const
         {
            return Base::operator[](vertex_id);
         }

         //! Access edge properties
         EdgeProps& operator[](const edge_descriptor& edge_id)
         {
#ifdef NDEBUG
            return edge_properties_[Base::operator[](edge_id)];
#else
            return edge_properties_.at(Base::operator[](edge_id));
#endif
         }

         //! Access edge properties
         const EdgeProps& operator[](
            const edge_descriptor& edge_id) const
         {
#ifdef NDEBUG
            return edge_properties_[Base::operator[](edge_id)];
#else
            return edge_properties_.at(Base::operator[](edge_id));
#endif
         }

         //! Directly access a graph bundle
         graph_bundled& operator[](boost::graph_bundle_t) {
            return Base::operator[](boost::graph_bundle);
         }

         const graph_bundled& operator[](boost::graph_bundle_t) const {
            return Base::operator[](boost::graph_bundle);
         }

         struct Potentials {
            using MagneticVectorPotential =
               typename Analysis::MagneticVectorPotential;
            using Iter = boost::transform_iterator<
               Potentials, boost::counting_iterator<vertex_descriptor> >;

            Potentials(const GridGraph& grid_graph)
               : grid_graph_(grid_graph)
            {
            }

            MagneticVectorPotential operator()(const vertex_descriptor& v) const
            {
               return grid_graph_[v]->potential();
            };

         private:
            const GridGraph& grid_graph_;
         };

         struct Positions {
            using Point = CartesianSystem::Point;
            using Iter = boost::transform_iterator<
               Positions, boost::counting_iterator<vertex_descriptor> >;

            Positions(const GridGraph& grid_graph)
               : grid_graph_(grid_graph)
            {
            }

            Point operator()(const vertex_descriptor& v) const
            {
               return grid_graph_[v]->position();
            }

         private:
            const GridGraph& grid_graph_;
         };

         template <typename Item, typename Iter = typename Item::Iter>
         Iter begin() const {
            return boost::make_transform_iterator(
               vertices(*this).first, Item(*this));
         }

         template <typename Item, typename Iter = typename Item::Iter>
         Iter end() const {
            return boost::make_transform_iterator(
               vertices(*this).second, Item(*this));
         }

      private:
         EdgePropertyVector edge_properties_;
      };
   }
}

//--------------------------------------------------------------
//
//   Function implementations
//
//--------------------------------------------------------------

#include <numeric>

#include "vertex.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Check that a GridGraphPointer points to something
      /*!
       * \tparam Analysis is the analysis type
       * \param gpointer is the pointer to check
       * \throw InvalidArgument in the case the pointer is nullptr
       * \relates GridGraph
       */
      template<typename GridGraphPointer>
      inline void CheckValidGridGraphPointer(GridGraphPointer ggp)
      {
#ifndef NDEBUG
         if (nullptr == ggp)
            throw InvalidArgument("Invalid GridGraphPointer!");
#endif
      }
   }

   //! Retrieve number of equations
   /*!
    * @param ggp points to a GridGraph
    *
    * @return number of equations
    *
    * \relates GridGraph
    */
   template<typename GridGraphPointer>
   inline EquationNumber NumberOfEquations(GridGraphPointer ggp)
   {
      return ((*ggp)[boost::graph_bundle_t()]).number_of_equations();
   }

   //! Set all vertices' potential from the result vector
   /*!
    * @param ggp points to a GridGraph
    *
    * @param potentials is the equation system's result vector
    *
    * \relates GridGraph
    */
   template<typename GridGraphPointer,
            typename Analysis =
            typename GridGraphPointer::element_type::AnalysisT>
   inline void SetPotential(
      GridGraphPointer ggp,
      const typename Analysis::RhsVector& potentials)
   {
      for (auto v = vertices(*ggp); v.first != v.second; ++(v.first)) {
         (*ggp)[*(v.first)]->set_potential(potentials);
      }
   }

   //! Set the flux density for all edges in a grid graph
   /*!
    * Call this function after the potentials have been set.
    *
    * @param ggp points to a GridGraph
    *
    * \bug Find some clever algorithm to determine the flux density's
    * sign.
    *
    * \relates GridGraph
    */
   template <typename GridGraphPointer>
   void UpdateMagneticFluxDensity(GridGraphPointer ggp)
   {
      for (auto e = edges(*ggp); e.first != e.second; ++(e.first)) {
         const auto source_vertex = source(*(e.first), *ggp);
         const auto target_vertex = target(*(e.first), *ggp);
         if (source_vertex < target_vertex) {
            const auto a1 = (*ggp)[source_vertex]->potential();
            const auto a2 = (*ggp)[target_vertex]->potential();
            const Length edge_length = (*ggp)[*(e.first)]->length();
            (*ggp)[*(e.first)]->set_flux_density((a2 - a1) / edge_length);
         }
      }
   }

   //! Switch to the next time step
   /*!
    *  \relates GridGraph
    */
   template <typename GridGraphPointer>
   void NextTimeStep(GridGraphPointer grid_graph)
   {
      (*grid_graph)[boost::graph_bundle_t()].NextTimeStep();
   }

   //! Find the vertex nearest to a given point
   /*!
    * @param ggp points to a GridGraph
    *
    * @param p is the approximate vertex position
    *
    * @returns the vertex closest to point
    *
    * @relates GridGraph
    */
   template <typename GridGraphPointer>
   typename GridGraphPointer::element_type::vertex_descriptor
   VertexNearPoint(GridGraphPointer ggp,
                   const CartesianSystem::Point& p)
   {
      using vertex_descriptor =
         typename GridGraphPointer::element_type::vertex_descriptor;

      vertex_descriptor current_vertex = 0, new_vertex=0;
      Length current_distance = CartesianSystem::EuclideanDistance(
         ((*ggp)[current_vertex])->position(), p);

      do {
         current_vertex = new_vertex;

         for (auto edges = out_edges(current_vertex, *ggp);
              edges.first != edges.second; ++edges.first) {
            const vertex_descriptor tv
               = target(*(edges.first), *ggp);
            const Length new_distance = CartesianSystem::EuclideanDistance(
               ((*ggp)[tv])->position(), p);
            if (new_distance < current_distance) {
               new_vertex = tv;
               current_distance = new_distance;
            }
         }
      } while (current_vertex != new_vertex);

      return new_vertex;
   }

   //! Write out as graphviz
   /*!
    * @param ggp points to a grid graph
    *
    * @param filename is the file to write to
    *
    * @relates GridGraph
    */
   template <typename GridGraphPointer>
   void WriteToGraphviz(GridGraphPointer ggp, const char* filename)
   {
      std::ofstream graphviz_file(filename);

      boost::write_graphviz(graphviz_file, *ggp);
   }

   //! A path inside the grid graph
   /*!
    *  \relates GridGraph
    */
   template <typename GridGraphPointer>
   struct PathInGridGraph
   {
   protected:
      using GridGraphType = typename GridGraphPointer::element_type;
      using Analysis = typename GridGraphType::AnalysisT;
      using MagneticFluxDensity = typename Analysis::MagneticFluxDensity;
      using MagneticVectorPotential =
         typename Analysis::MagneticVectorPotential;
      using MagneticVoltage = typename Analysis::MagneticVoltage;
      using VertexDescriptor =
         typename boost::graph_traits<GridGraphType>::vertex_descriptor;
      using EdgeDescriptor =
         typename boost::graph_traits<GridGraphType>::edge_descriptor;
      using Point = CartesianSystem::Point;
      using VertexVector = std::vector<VertexDescriptor>;
      using EdgeVector = std::vector<EdgeDescriptor>;
      using VectorNormalized = CartesianSystem::VectorNormalized;

   public:
      //! Initialization
      /*!
       *  \param grid_graph pointer to a grid graph
       *
       *  \param start desired starting point.  The real starting
       *  point is a vertex which is as close as possible to start.
       *
       *  \param end desired end point.  The real end
       *  point is a vertex which is as close as possible to end.
       *
       *  \throws RuntimeError if the resulting path definition is
       *  invalid.
       */
      template <typename PointIter>
      PathInGridGraph(GridGraphPointer grid_graph,
                      PointIter start, PointIter end)
         : grid_graph_(grid_graph)
      {
         vertices_.emplace_back(VertexNearPoint(grid_graph, *start));
         while(++start != end)
            this->AddSection(*start);
         vertices_.shrink_to_fit();
         this->CheckValid();
      }

      //! Determine the average flux density along this path
      MagneticFluxDensity average_flux_density() const;

      //! Determine the magnetic voltage
      MagneticVoltage magnetic_voltage() const;

      //! Total length
      Length total_length() const;

      struct Potentials {
         using Iter = boost::transform_iterator<
            Potentials, typename VertexVector::const_iterator>;

         Potentials(GridGraphPointer grid_graph)
            : grid_graph_(grid_graph)
         {
         }

         MagneticVectorPotential operator()(const VertexDescriptor& v) const
         {
            return (*(this->grid_graph_))[v]->potential();
         };

      private:
         GridGraphPointer grid_graph_;
      };

      struct Positions {
         using Iter = boost::transform_iterator<
            Positions, typename VertexVector::const_iterator>;

         Positions(GridGraphPointer grid_graph)
            : grid_graph_(grid_graph)
         {
         }

         Point operator()(const VertexDescriptor& v) const
         {
            return (*(this->grid_graph_))[v]->position();
         }

      private:
         GridGraphPointer grid_graph_;
      };

      template <typename Item, typename Iter = typename Item::Iter>
      Iter begin() const {
         return boost::make_transform_iterator(
            vertices_.begin(), Item(grid_graph_));
      }

      template <typename Item, typename Iter = typename Item::Iter>
      Iter end() const {
         return boost::make_transform_iterator(
            vertices_.end(), Item(grid_graph_));
      }

   protected:
      //! Do the initialization
      /*!
       *  Find a path from start to end.
       */
      void AddSection(const Point&);

      //! Check that the final path definition is valid
      void CheckValid() const
      {
         if (vertices_.size() < 2 or vertices_[0] == vertices_[1])
            throw RuntimeError("Invalid path definition");
      }

   private:
      VertexVector vertices_;
      EdgeVector edges_;
      GridGraphPointer grid_graph_;
   };

   template <typename GridGraphPointer>
   inline void PathInGridGraph<GridGraphPointer>::AddSection(
      const Point& end_point)
   {
      VertexDescriptor current_vertex=vertices_.back();
      VertexDescriptor new_vertex=current_vertex;
      Length current_distance = CartesianSystem::EuclideanDistance(
         ((*grid_graph_)[current_vertex])->position(), end_point);
      bool done = false;

      do {
         current_vertex = new_vertex;

         EdgeDescriptor current_edge;

         for (auto edges = out_edges(current_vertex, *grid_graph_);
              edges.first != edges.second; ++edges.first) {
            const VertexDescriptor target_vertex
               = target(*(edges.first), *grid_graph_);
            const Length new_distance = CartesianSystem::EuclideanDistance(
               ((*grid_graph_)[target_vertex])->position(), end_point);
            if (new_distance < current_distance) {
               new_vertex = target_vertex;
               current_distance = new_distance;
               current_edge = *(edges.first);
            }
         }

         done = (current_vertex == new_vertex);

         if (not done) {
            vertices_.emplace_back(new_vertex);
            edges_.emplace_back(current_edge);
         }
      } while (not done);
   }

   template <typename GridGraphPointer>
   inline Length PathInGridGraph<GridGraphPointer>::total_length() const
   {
      using boost::units::si::metre;

      Length total_length = .0*metre;

      for (auto iter=edges_.begin(), end=edges_.end(); iter != end; ++iter)
         total_length += (*grid_graph_)[*iter]->length();

      return total_length;
   }

   template <typename GridGraphPointer>
   inline auto
   PathInGridGraph<GridGraphPointer>::average_flux_density() const
      -> MagneticFluxDensity
   {
      return ((*grid_graph_)[vertices_.front()]->potential()
              -(*grid_graph_)[vertices_.back()]->potential())/total_length();
   }

   template <typename GridGraphPointer>
   inline auto
   PathInGridGraph<GridGraphPointer>::magnetic_voltage() const
      -> MagneticVoltage
   {
      // using boost::units::si::square_metre;

      MagneticVoltage voltage_sum = MagneticVoltage::from_value(0.);

      if (vertices_.size() < 3) return voltage_sum;

      for (auto iter = (vertices_.begin())+1, end = (vertices_.end())-1;
           iter != end; ++iter) {

         unsigned short dangling_edge_count = 0;
         MagneticVoltage this_vertex_voltage_sum
            = MagneticVoltage::from_value(0.);

         for (auto all_edges=out_edges(*iter, *grid_graph_);
              all_edges.first != all_edges.second; ++(all_edges.first)) {

            auto target_vertex = target(*(all_edges.first), *grid_graph_);

            // ignore edges onto the path
            if (target_vertex == *(iter-1) or
                target_vertex == *(iter+1)) continue;

            // determine whether vertex is left or right from path
            const auto left_or_right_from_path =
               Signum(
                  VectorProduct(
                     Substract(
                        (*grid_graph_)[target_vertex]->position(),
                        (*grid_graph_)[*iter]->position()),
                     Substract(
                        (*grid_graph_)[*(iter+1)]->position(),
                        (*grid_graph_)[*iter]->position())));

            const MagneticVoltage helper1 =
               ((*grid_graph_)[*iter]->potential()
                -(*grid_graph_)[target_vertex]->potential())
               *(*grid_graph_)[*(all_edges.first)]->reluctivity();

            // Sum up magnetic voltage of all "dangling" edges
            this_vertex_voltage_sum +=  MagneticVoltage::from_value(
               helper1.value()
               *(*grid_graph_)[*(all_edges.first)]->length_coefficient()
               *(.0+left_or_right_from_path));

            ++dangling_edge_count;
         }

         voltage_sum +=
            MagneticVoltage::from_value(
               this_vertex_voltage_sum.value()/(.0+dangling_edge_count));
      }

      return voltage_sum;
   }

#ifndef INSTANTIATE_
   // extern template void UpdateMagneticFluxDensity(GridGraphPointer<TimeTransientAnalysis<> > );
   // extern template void UpdateMagneticFluxDensity(GridGraphPointer<TimeHarmonicAnalysis<> > );
#endif
}

#endif // GRID_GRAPH_HPP_YV9h1g37_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
