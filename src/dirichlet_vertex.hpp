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

/*!
 *  \file
 *
 *  \brief Implement the dirichlet boundary condition.
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef DIRICHLET_VERTEX_HPP_nrSrK82d_
#define DIRICHLET_VERTEX_HPP_nrSrK82d_

#include "config.h"

#include "grid_graph-fwd.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Vertex used for modeling the Dirichlet boundary condition
      /*!
       *  This vertex type is used to model the Dirichlet boundary
       *  condition.  This boundary condition forces the vertex'
       *  magnetic vector potential to a fixed value.
       *
       *  @tparam Analysis defines the analysis type (time transient
       *  or time harmonic).
       */
      template <typename Analysis, template <typename> class Traits>
      struct DirichletVertex :
         Vertex<Analysis, Traits>
      {
      protected:
         using Base = Vertex<Analysis, Traits>;
         using VertexTraits = Traits<Analysis>;
         using StandardVertexType = StandardVertex<Analysis, Traits>;
         using NeumannVertexType = NeumannVertex<Analysis, Traits>;
         using GridGraphPtr = typename Base::GridGraphPtr;
         using AssemblingMatrix = typename Base::AssemblingMatrix;
         using CoeffCoeffDerivedPair = typename Base::CoeffCoeffDerivedPair;
         using CoefficientQuantity = typename Base::CoefficientQuantity;
         using MagneticVectorPotential = typename Base::MagneticVectorPotential;
         using Point = typename Base::Point;
         using RhsVector = typename Base::RhsVector;
         using VertexDescriptor = typename Base::VertexDescriptor;
         using Time = typename Analysis::Time;
         using EdgeProps = typename Base::EdgeProps;
         using PotentialFunction =
            std::function<MagneticVectorPotential (const Point&, Time)>;

      public:
         //! Initialize a vertex for a homogeneous dirichlet condition
         DirichletVertex()
            : Base(), current_time_(Time::from_value(0.)),
              potential_function_(
                 [](const Point&, Time) -> MagneticVectorPotential {
                    return MagneticVectorPotential::from_value(.0);
                 })
         {
         }

         //! Initialize a vertex for a homogeneous dirichlet condition
         /*!
          * @param position is the vertex's position within the global
          * cartesian coordinate system
          */
         explicit DirichletVertex(const Point& position)
            : Base(position), current_time_(Time::from_value(.0)),
              potential_function_(
                 [](const Point&, Time) -> MagneticVectorPotential {
                    return MagneticVectorPotential::from_value(.0);
                 })
         {
         }

         //! Initialization
         DirichletVertex(const Point& position,
                         PassType<MagneticVectorPotential> fixed_potential)
            : Base(position), current_time_(Time::from_value(.0)),
              potential_function_(
                 [fixed_potential](const Point&, Time) -> MagneticVectorPotential {
                    return fixed_potential;
                 })
         {
         }

         //! Initialization
         /*!
          * @param position is the vertex' position (global coordinate
          * system).
          *
          * @param potential_function takes a Point (cartesian) and a
          * time/frequency parameter and computes the related
          * potential.
          */
         DirichletVertex(
            const Point& position,
            const PotentialFunction& potential_function)
            : Base(position), current_time_(Time::from_value(.0)),
              potential_function_(potential_function)
         {
         }

         //! Destruction
         virtual ~DirichletVertex()
         {
         }

         void set_potential(const RhsVector&) override final
         {
         }

         //! Set current time
         void set_current_time(Time t) const
         {
            current_time_ = t;
         }

         void AddToEquationSystem(
            VertexDescriptor vertex_id, GridGraphPtr grid_graph,
            AssemblingMatrix& system_triplets, RhsVector& rhs) const override
         {
            using boost::out_edges;

            current_time_ = (*grid_graph)[boost::graph_bundle].current_time();

            for (auto e = out_edges(vertex_id, *grid_graph);
                 e.first != e.second; ++(e.first))
            {
               (*grid_graph)[target(*(e.first), *grid_graph)]->AddEdge(
                  this, (*grid_graph)[*(e.first)], system_triplets, rhs);
            }
         }

         void AddToEquationSystem(
            VertexDescriptor vertex_id, GridGraphPtr grid_graph,
            AssemblingMatrix& system_triplets, AssemblingMatrix&,
            RhsVector& rhs) const override
         {
            this->AddToEquationSystem(
               vertex_id, grid_graph, system_triplets, rhs);
         }

      protected:
         CoefficientQuantity AddEdge(
            const StandardVertexType*, const EdgeProps& edge_properties,
            AssemblingMatrix&, RhsVector&) const override
         {
            return edge_properties->SystemMatrixCoefficient();
         }

         CoeffCoeffDerivedPair AddEdge(
            const StandardVertexType*, const EdgeProps& edge_properties,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const override
         {
            return edge_properties->SystemAndJacobiMatrixCoefficient();
         }

         CoefficientQuantity AddEdge(
            const NeumannVertexType*, const EdgeProps& edge_properties,
            AssemblingMatrix&, RhsVector&) const override
         {
            return edge_properties->SystemMatrixCoefficient();
         }

         CoeffCoeffDerivedPair AddEdge(
            const NeumannVertexType*, const EdgeProps& edge_properties,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const override
         {
            return edge_properties->SystemAndJacobiMatrixCoefficient();
         }

         CoefficientQuantity AddEdge(
            const DirichletVertex*, const EdgeProps&,
            AssemblingMatrix&, RhsVector&) const override
         {
            return CoefficientQuantity::from_value(0);
         }

         CoeffCoeffDerivedPair AddEdge(
            const DirichletVertex*, const EdgeProps&,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const override
         {
            return std::make_tuple(CoefficientQuantity::from_value(0),
                                   CoefficientQuantity::from_value(0));
         }

      protected:
         MagneticVectorPotential DoPotential() const override
         {
            return potential_function_(this->position(), current_time_);
         }

      private:
         mutable Time current_time_;
         PotentialFunction potential_function_;
      };
   }
}

#endif // DIRICHLET_VERTEX_HPP_nrSrK82d_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
