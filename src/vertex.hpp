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
 *  \brief Vertex base class
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef VERTEX_HPP_iteW0Wv5_
#define VERTEX_HPP_iteW0Wv5_

#include <memory>

#include "config.h"

#ifdef HAVE_LIBLOKI
#include <loki/SmallObj.h>
#endif

#include "boost_units_extra.hpp"
#include "coordinate_systems.hpp"
#include "types.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Vertex base class
      /*!
       * \tparam Analysis is the analysis type
       */
      template<typename Analysis, template <typename> class Traits>
      struct Vertex
#ifdef HAVE_LIBLOKI
         : Loki::SmallObject<>
#endif
      {
      protected:
         //! Standard vertex type
         using StandardVertexType = StandardVertex<Analysis, Traits>;

         //! Dirichlet vertex type
         using DirichletVertexType = DirichletVertex<Analysis, Traits>;

         //! Neumann vertex type
         using NeumannVertexType = NeumannVertex<Analysis, Traits>;

         using VertexTraits = Traits<Analysis>;

         //! Grid graph pointer
         using GridGraphPtr = typename VertexTraits::GridGraphPtr;

         //! System matrix coefficient type
         using CoefficientQuantity = typename VertexTraits::CoefficientQuantity;

         //! The equation system's right hand side
         using RhsVector = typename Analysis::RhsVector;

         //! Vertex descriptor type
         using VertexDescriptor = typename VertexTraits::VertexDescriptor;

         //! Magnetic vector potential type
         using MagneticVectorPotential =
            typename Analysis::MagneticVectorPotential;

         //! EdgeProperties
         using EdgeProps = typename VertexTraits::EdgeProps;

         //! Assembling matrix type
         using AssemblingMatrix = typename Analysis::AssemblingMatrix;

         //! Finite difference and jacobi coefficient
         using CoeffCoeffDerivedPair =
            std::tuple<CoefficientQuantity, CoefficientQuantity>;

         //! Point type
         using Point = CartesianSystem::Point;

         //! Type for passing a magnetic vector potential value
         using MagneticVectorPotentialPassType =
            typename detail::PassType<MagneticVectorPotential>;

      public:
         //! Default initialization
         Vertex()
            : position_(
               std::make_tuple(.0 * boost::units::si::metre,
                               .0 * boost::units::si::metre))
         {
         }

         //! Initialization
         /*!
          * @param position is the vertex's position within the global
          * cartesian coordinate system
          */
         explicit Vertex(const Point& position) :
            position_(position)
         {
         }

         //! Destruction
         virtual ~Vertex()
         {
         }

         //! Get the vertex's position
         /*!
          * @return the vertex's position in the global coordinate system
          */
         const Point& position() const
         {
            return position_;
         }

         //! Set the vertex's position
         /*!
          * @param position is the new position
          */
         void set_position(const Point& position)
         {
            position_ = position;
         }

         //! Set the vertex's potential from the result vector.
         /*!
          * @param rhs is the equation system's result vector
          */
         virtual void set_potential(const RhsVector&) =0;

         //! Retrieve the vertex' potential
         MagneticVectorPotential potential() const {
            return this->DoPotential();
         }

      protected:
         //! Do retrieve the vertex' potential
         virtual MagneticVectorPotential DoPotential() const =0;

      public:
         //! Add this vertex to a linear equation system
         /*!
          * This is stage one in a double dispatch process.
          *
          * @param vertex represents this vertex within the grid graph.
          *
          * @param grid_graph is the grid graph.
          *
          * @param system_triplets is a vector for collecting all system
          * matrix entries.
          *
          * @param rhs is the equation system's right hand side (current
          * density vector).
          */
         virtual void AddToEquationSystem(
            VertexDescriptor, GridGraphPtr,
            AssemblingMatrix&, RhsVector&) const =0;

         //! Add this vertex to a nonlinear equation system
         /*!
          * This is stage one in a double dispatch process.
          *
          * @param vertex represents this vertex within the grid graph.
          * @param grid_graph is the grid graph.
          *
          * @param system_triplets is a vector for collecting all system
          * matrix entries.
          *
          * @param jacobi_triplets is a vector for collecting all jacobi
          * matrix entries.
          *
          * @param rhs is the equation system's right hand side (current
          * density vector).
          */
         virtual void AddToEquationSystem(
            VertexDescriptor, GridGraphPtr,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const =0;

         //! Add a single edge to a linear equation system.
         /*!
          * This is stage two in a double dispatch process.
          *
          * @param source is the edge's source vertex (type)
          * @param edge is the edge which shall be added.
          * @param grid_graph describes the investigated problem
          *
          * @param system_triplets is a vector for collecting all system
          * matrix entries.
          *
          * @param rhs is the equation system's right hand side (current
          * density vector).
          *
          * @return the inserted coefficient
          */
         virtual CoefficientQuantity AddEdge(
            const StandardVertexType*, const EdgeProps&,
            AssemblingMatrix&, RhsVector&) const =0;

         //! Add a single edge to a nonlinear equation system.
         /*!
          * This is stage two in a double dispatch process.
          *
          * @param source is the edge's source vertex (type)
          * @param edge is the edge which shall be added.
          * @param grid_graph describes the investigated problem
          *
          * @param system_triplets is a vector for collecting all system
          * matrix entries.
          *
          * @param jacobi_triplets is a vector for collecting all jacobi
          * matrix entries.
          *
          * @param rhs is the equation system's right hand side (current
          * density vector).
          *
          * @return the inserted coefficient
          */
         virtual CoeffCoeffDerivedPair AddEdge(
            const StandardVertexType*, const EdgeProps&,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const =0;

         //! Add a single edge to a linear equation system.
         /*!
          * This is stage two in a double dispatch process.
          *
          * @param source is the edge's source vertex (type)
          * @param edge is the edge which shall be added.
          * @param grid_graph describes the investigated problem
          *
          * @param system_triplets is a vector for collecting all system
          * matrix entries.
          *
          * @param rhs is the equation system's right hand side (current
          * density vector).
          *
          * @return the inserted coefficient
          */
         virtual CoefficientQuantity AddEdge(
            const NeumannVertexType*, const EdgeProps&,
            AssemblingMatrix&, RhsVector&) const =0;

         //! Add a single edge to a nonlinear equation system.
         /*!
          * This is stage two in a double dispatch process.
          *
          * @param source is the edge's source vertex (type)
          * @param edge is the edge which shall be added.
          * @param grid_graph describes the investigated problem
          *
          * @param system_triplets is a vector for collecting all system
          * matrix entries.
          *
          * @param jacobi_triplets is a vector for collecting all jacobi
          * matrix entries.
          *
          * @param rhs is the equation system's right hand side (current
          * density vector).
          *
          * @return the inserted coefficient
          */
         virtual CoeffCoeffDerivedPair AddEdge(
            const NeumannVertexType*, const EdgeProps&,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const =0;

         //! Add a single edge to a linear equation system.
         /*!
          * This is stage two in a double dispatch process.
          *
          * @param source is the edge's source vertex (type)
          *
          * @param edge is the edge which shall be added.
          *
          * @param grid_graph describes the investigated problem
          *
          * @param system_triplets is a vector for collecting all system
          * matrix entries.
          *
          * @param rhs is the equation system's right hand side (current
          * density vector).
          *
          * @return the inserted coefficient
          */
         virtual CoefficientQuantity AddEdge(
            const DirichletVertexType*, const EdgeProps&,
            AssemblingMatrix&, RhsVector&) const =0;

         //! Add a single edge to a nonlinear equation system.
         /*!
          * This is stage two in a double dispatch process.
          *
          * @param source is the edge's source vertex (type)
          * @param edge is the edge which shall be added.
          * @param grid_graph describes the investigated problem
          *
          * @param system_triplets is a vector for collecting all system
          * matrix entries.
          *
          * @param jacobi_triplets is a vector for collecting all jacobi
          * matrix entries.
          *
          * @param rhs is the equation system's right hand side (current
          * density vector).
          *
          * @return the inserted coefficient
          */
         virtual CoeffCoeffDerivedPair AddEdge(
            const DirichletVertexType*, const EdgeProps&,
            AssemblingMatrix&, AssemblingMatrix&, RhsVector&) const =0;

      private:
         Point position_;
      };
   }
}

#endif // VERTEX_HPP_iteW0Wv5_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
