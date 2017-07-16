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
 *  \brief Standard vertex (inner vertex and homogeneous neumann condition)
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef STANDARD_VERTEX_HPP_SaFtkO8R_
#define STANDARD_VERTEX_HPP_SaFtkO8R_

#include "config.h"

#include <boost/array.hpp>

#include "analysis.hpp"
#include "grid_graph-fwd.hpp"
#include "vertex_properties.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! General purpose vertex
      /*!
       *  \tparam Analysis is the analysis type
       */
      template<typename Analysis, template <typename> class Traits>
      struct StandardVertex :
         Vertex<Analysis, Traits>,
         Traits<Analysis>::StandardVertexProperties
      {
      protected:
         //! Base class
         using Base = Vertex<Analysis, Traits>;
         using ThisType = StandardVertex<Analysis, Traits>;

         using NeumannVertexType = NeumannVertex<Analysis, Traits>;
         using DirichletVertexType = DirichletVertex<Analysis, Traits>;
         using VertexTraits = Traits<Analysis>;
         // using GridGraph = typename VertexTraits::GridGraph;
         using GridGraphPtr = typename Base::GridGraphPtr;
         using AssemblingMatrix = typename Base::AssemblingMatrix;
         using CoeffCoeffDerivedPair = typename Base::CoeffCoeffDerivedPair;
         using CoefficientQuantity = typename Base::CoefficientQuantity;
         using MagneticVectorPotential = typename Base::MagneticVectorPotential;
         using Point = typename Base::Point;
         using RhsVector = typename Base::RhsVector;
         using VertexDescriptor = typename Base::VertexDescriptor;
         using StandardProperties =
            typename VertexTraits::StandardVertexProperties;
         using EdgeProps = typename Base::EdgeProps;

      public:
         //! Default initialization
         StandardVertex() : Base(), StandardProperties(), equation_number_(0)
         {
         }

         //! Initialization
         /*!
          * @param position is the vertex' position within the global
          * cartesian coordinate system
          *
          * @param arguments are forwarded to the properties object
          * this class is derived from.
          */
         template <typename ... Args>
         StandardVertex(const Point& position, Args ... arguments)
            : Base(position), StandardProperties(arguments ...),
              equation_number_(0)
         {
         }

         //! Destruction
         virtual ~StandardVertex()
         {
         }

         //! Get the vertex's equation number
         /*!
          * @return the vertex's equation number
          */
         EquationNumber equation_number() const
         {
            return equation_number_;
         }

         //! Set the vertex's equation number
         /*!
          * @param equationNumber is the new equation number
          */
         void set_equation_number(EquationNumber equation_number)
         {
            equation_number_ = equation_number;
         }

         MagneticVectorPotential DoPotential() const override
         {
            return this->DoGetPotential();
         }

         // void set_potential(MagneticVectorPotential new_potential)
         // {
         //    potentials_.SetPotential(new_potential);
         // }

         void set_potential(const RhsVector& rhs) override final
         {
            using boost::units::si::weber_per_metre;

            this->DoSetPotential(
               rhs[this->equation_number_] * weber_per_metre);
         }

         void AddToEquationSystem(
            VertexDescriptor, GridGraphPtr,
            AssemblingMatrix&, RhsVector&) const override;
         void AddToEquationSystem(
            VertexDescriptor, GridGraphPtr,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const override;

      protected:
         CoefficientQuantity AddEdge(
            const StandardVertex*, const EdgeProps&,
            AssemblingMatrix&, RhsVector&) const override;
         CoeffCoeffDerivedPair AddEdge(
            const StandardVertex*, const EdgeProps&,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const override;
         CoefficientQuantity AddEdge(
            const NeumannVertexType*, const EdgeProps&,
            AssemblingMatrix&, RhsVector&) const override;
         CoeffCoeffDerivedPair AddEdge(
            const NeumannVertexType*, const EdgeProps&,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const override;
         CoefficientQuantity AddEdge(
            const DirichletVertexType*, const EdgeProps&,
            AssemblingMatrix&, RhsVector&) const override;
         CoeffCoeffDerivedPair AddEdge(
            const DirichletVertexType*, const EdgeProps&,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const override;

      private:
         EquationNumber equation_number_;
      };
   }
}

//--------------------------------------------------------------
//
//   Function implementations
//
//--------------------------------------------------------------

#include "dirichlet_vertex.hpp"
#include "neumann_vertex.hpp"

namespace PACKAGE_NS {
   namespace detail {
      template <typename Analysis, template <typename> class Traits>
      auto StandardVertex<Analysis, Traits>::AddEdge(
         const DirichletVertexType* source,
         const EdgeProps& edge_properties,
         AssemblingMatrix& system_triplets, RhsVector& rhs) const
         -> CoefficientQuantity
      {
         // using DataType = typename Analysis::DataType;

         const auto coefficient = edge_properties->SystemMatrixCoefficient();

         rhs(this->equation_number()) += in_base_units(
            (coefficient * source->potential()));

         return coefficient;
      }

      template <typename Analysis, template <typename> class Traits>
      auto StandardVertex<Analysis, Traits>::AddEdge(
         const DirichletVertexType* source,
         const EdgeProps& edge_properties,
         AssemblingMatrix& system_triplets, AssemblingMatrix& jacobi_triplets,
         RhsVector& rhs) const -> CoeffCoeffDerivedPair
      {
         // using DataType = typename Analysis::DataType;

         const auto fd_jacobi_coefficient
            = edge_properties->SystemAndJacobiMatrixCoefficient();

         rhs(this->equation_number()) += in_base_units(
            (std::get<0>(fd_jacobi_coefficient) * source->potential()));

         return fd_jacobi_coefficient;
      }

      template <typename Analysis, template <typename> class Traits>
      auto StandardVertex<Analysis, Traits>::AddEdge(
         const StandardVertex* source,
         const EdgeProps& edge_properties,
         AssemblingMatrix& system_triplets, RhsVector& rhs) const
         -> CoefficientQuantity
      {
         const EquationNumber
            source_equation_number = source->equation_number(),
            target_equation_number = this->equation_number();

         const auto fd_coefficient = edge_properties->SystemMatrixCoefficient();

         system_triplets.emplace_back(
            source_equation_number, target_equation_number,
            in_base_units(-fd_coefficient));

         return fd_coefficient;
      }

      template <typename Analysis, template <typename> class Traits>
      auto StandardVertex<Analysis, Traits>::AddEdge(
         const StandardVertex* source,
         const EdgeProps& edge_properties,
         AssemblingMatrix& system_triplets, AssemblingMatrix& jacobi_triplets,
         RhsVector& rhs) const -> CoeffCoeffDerivedPair
      {
         const EquationNumber
            source_equation_number = source->equation_number(),
            target_equation_number = this->equation_number();

         const auto fd_jacobi_coefficient
            = edge_properties->SystemAndJacobiMatrixCoefficient();

         system_triplets.emplace_back(
            source_equation_number, target_equation_number,
            in_base_units(-std::get<0>(fd_jacobi_coefficient)));
         jacobi_triplets.emplace_back(
            source_equation_number, target_equation_number,
            in_base_units(-std::get<1>(fd_jacobi_coefficient)));

         return fd_jacobi_coefficient;
      }

      template <typename Analysis, template <typename> class Traits>
      auto StandardVertex<Analysis, Traits>::AddEdge(
         const NeumannVertexType* source,
         const EdgeProps& edge_properties,
         AssemblingMatrix& system_triplets, RhsVector& rhs) const
         -> CoefficientQuantity
      {
         return this->AddEdge(
            static_cast<const ThisType*>(source), edge_properties,
            system_triplets, rhs);
      }

      template <typename Analysis, template <typename> class Traits>
      auto StandardVertex<Analysis, Traits>::AddEdge(
         const NeumannVertexType* source,
         const EdgeProps& edge_properties,
         AssemblingMatrix& system_triplets, AssemblingMatrix& jacobi_triplets,
         RhsVector& rhs) const -> CoeffCoeffDerivedPair
      {
         return this->AddEdge(
            static_cast<const ThisType*>(source), edge_properties,
            system_triplets, jacobi_triplets, rhs);
      }

      template <typename Analysis, template <typename> class Traits>
      void StandardVertex<Analysis, Traits>::AddToEquationSystem(
         VertexDescriptor vertex_id, GridGraphPtr grid_graph,
         AssemblingMatrix& system_triplets, RhsVector& rhs) const
      {
         using DataType = typename Analysis::DataType;
         using Time = typename Analysis::Time;

         const Time
            current_time = (*grid_graph)[boost::graph_bundle].current_time();

         rhs(this->equation_number_) +=
            in_base_units(this->StandardRhs(current_time));

         DataType coefficient_sum = .0;

         for (auto e = out_edges(vertex_id, *grid_graph); e.first != e.second;
              ++(e.first)) {
            coefficient_sum += in_base_units(
               (*grid_graph)[target(*(e.first), *grid_graph)]->AddEdge(
                  this, (*grid_graph)[*(e.first)], system_triplets, rhs));
         }

         system_triplets.emplace_back(
            this->equation_number_, this->equation_number_, coefficient_sum);
      }

      template <typename Analysis, template <typename> class Traits>
      void StandardVertex<Analysis, Traits>::AddToEquationSystem(
         VertexDescriptor vertex_id, GridGraphPtr grid_graph,
         AssemblingMatrix& system_triplets, AssemblingMatrix& jacobi_triplets,
         RhsVector& rhs) const
      {
         using DataType = typename Analysis::DataType;
         using std::get;
         using Time = typename Analysis::Time;

         const Time
            current_time = (*grid_graph)[boost::graph_bundle].current_time();

         rhs(this->equation_number_) +=
            in_base_units(this->StandardRhs(current_time));

         DataType coefficient_sum = .0;
         DataType jacobi_coefficient_sum = .0;

         for (auto e = out_edges(vertex_id, *grid_graph); e.first != e.second;
              ++(e.first))
         {
            const auto fd_jacobi_coefficient
               = (*grid_graph)[target(*(e.first), *grid_graph)]->AddEdge(
                  this, (*grid_graph)[*(e.first)],
                  system_triplets, jacobi_triplets, rhs);

            coefficient_sum += in_base_units(get<0>(fd_jacobi_coefficient));
            jacobi_coefficient_sum += in_base_units(
               get<1>(fd_jacobi_coefficient));
         }

         system_triplets.emplace_back(
            this->equation_number_, this->equation_number_, coefficient_sum);
         jacobi_triplets.emplace_back(
            this->equation_number_, this->equation_number_,
            jacobi_coefficient_sum);
      }
   }
}

#endif // STANDARD_VERTEX_HPP_SaFtkO8R_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
