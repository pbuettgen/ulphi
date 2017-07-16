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
 *  \brief Create vertex objects for vertices onto boundaries
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef BOUNDARY_VERTEX_FACTORY_HPP_tYB3hesK_
#define BOUNDARY_VERTEX_FACTORY_HPP_tYB3hesK_

#include "config.h"

#include "boundary-fwd.hpp"
#include "coordinate_systems.hpp"
#include "dirichlet_condition.hpp"
#include "grid_graph-fwd.hpp"
#include "indices.hpp"
#include "load_definition-fwd.hpp"
#include "loki_visitor_extra.hpp"
#include "neumann_condition.hpp"

#define NEW_ON_BOUNDARY(BOUNDARY, INDEX)		\
   void operator()(					\
      const BOUNDARY ## Boundary<Analysis>& boundary,	\
      INDEX ## Index i)					\
   {							\
      boundary_ = &boundary;				\
      i_ = boundary.AsBoundaryIndex(i);			\
      boundary.boundary_condition()->Accept(*this);	\
   }

namespace PACKAGE_NS {
   namespace detail {
      template <typename Analysis, template <typename> class VertexPtrType>
      struct BoundaryVertexFactoryBase :
	 Loki::BaseVisitor,
	 Loki::Visitor<DirichletCondition<Analysis>, void, true >,
	 Loki::Visitor<NeumannCondition<Analysis>, void, true >,
	 Loki::Visitor<VoidCondition<Analysis>, void, true >
      {
      private:
	 using VertexPtr = VertexPtrType<Analysis>;
	 using LoadDefPtr = ConstLoadDefinitionPointer<Analysis>;

      public:
	 //! Initialization
	 BoundaryVertexFactoryBase()
	    : boundary_(nullptr), i_(0), new_vertex_object_(nullptr)
	 {
	 }

	 //! Call a new vertex object into life
	 /*!
	  * @param boundary the boundary
	  *
	  * @param i is the vertex' index onto the boundary
	  */
	 void operator()(
	    const Boundary<Analysis>& boundary, BoundaryIndex i)
	 {
	    boundary_ = &boundary;
	    i_ = i;
	    boundary.boundary_condition()->Accept(*this);
	 }

	 NEW_ON_BOUNDARY(North, Column)
	 NEW_ON_BOUNDARY(West, Row)
	 NEW_ON_BOUNDARY(South, Column)
	 NEW_ON_BOUNDARY(East, Row)

	 void set_new_vertex_object(VertexPtr& new_vertex)
	 {
	    std::swap(new_vertex_object_, new_vertex);
	 }

	 VertexPtr& new_vertex_object() {
	    return new_vertex_object_;
	 }

      protected:
	 //! Vertex position within the global coordinate system
	 CartesianSystem::Point position() {
	    return boundary_->VertexPositionGlobal(i_);
	 }

	 //! Retrieve the load definition
	 LoadDefPtr load_definition() {
	    return boundary_->GetLoad(i_);
	 }

	 const Boundary<Analysis>* boundary_;
	 BoundaryIndex i_;
	 VertexPtr new_vertex_object_;
      };
   }

   namespace fi {

      //! Generate appropriate vertex objects for boundary vertices
      /*!
       * Informations from the grid and the boundary are combinded in
       * order to create an appropriate vertex object.
       *
       * @tparam Analysis is the analysis type
       */
      template <typename Analysis>
      struct BoundaryVertexFactory :
	 ::PACKAGE_NS::detail::BoundaryVertexFactoryBase<Analysis, VertexPointer>
      {
      private:
	 using VertexPtr = VertexPointer<Analysis>;
	 using LoadDefPtr = ConstLoadDefinitionPointer<Analysis>;

      public:
	 //! Visit a VoidCondition
	 void Visit(const VoidCondition<Analysis>&) override
	 {
	    this->new_vertex_object_ =
	       this->boundary_->HasBoundaryCoupling() ?
	       VertexPtr(nullptr) :
	       std::make_unique<StandardVertex<Analysis> >(
		  this->position(), this->load_definition(),
		  this->area());
	 }

	 //! Visit a DirichletCondition
	 void Visit(
	    const DirichletCondition<Analysis>& dirichlet_condition) override
	 {
	    this->new_vertex_object_ =
	       this->boundary_->HasBoundaryCoupling() ?
	       VertexPtr(nullptr) :
	       std::make_unique<DirichletVertex<Analysis> >(
		  this->position(), dirichlet_condition.potential_function());
	 }

	 //! Visit a NeumannCondition
	 void Visit(
	    const NeumannCondition<Analysis>& neumann_condition) override
	 {
	    this->new_vertex_object_ =
	       this->boundary_->HasBoundaryCoupling() ?
	       VertexPtr(nullptr) :
	       std::make_unique<NeumannVertex<Analysis> >(
		  this->position(),
		  neumann_condition.reluctivity_outside(),
		  neumann_condition.flux_density_outside(),
		  this->boundary_->EdgeLengthPerpendicular(this->i_),
		  this->boundary_->parallelity(),
		  this->load_definition(), this->area());
	 }

      protected:
	 //! Secondary grid's cell size
	 Area area() {
	    return this->boundary_->CellSize(this->i_);
	 }
      };
   }
}

#ifdef NEW_ON_BOUNDARY
#undef NEW_ON_BOUNDARY
#endif

#endif // BOUNDARY_VERTEX_FACTORY_HPP_tYB3hesK_
