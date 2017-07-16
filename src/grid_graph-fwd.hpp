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
 *  \brief Forward declarations for grid graph
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef GRID_GRAPHFWD_HPP_mRjov6k1_
#define GRID_GRAPHFWD_HPP_mRjov6k1_

#include <memory>

#include "config.h"

#include <boost/graph/compressed_sparse_row_graph.hpp>

namespace PACKAGE_NS {

   //
   // --- Vertex types -----------------------------------------------
   //
#define DEFINE_VERTEX_TYPES                                             \
   template <typename Analysis>                                         \
   using Vertex =                                                       \
      ::PACKAGE_NS::detail::Vertex<Analysis, VertexTraits>;             \
   template <typename Analysis>                                         \
   using StandardVertex =                                               \
      ::PACKAGE_NS::detail::StandardVertex<Analysis, VertexTraits>;     \
   template <typename Analysis>                                         \
   using NeumannVertex =                                                \
      ::PACKAGE_NS::detail::NeumannVertex<Analysis, VertexTraits>;      \
   template <typename Analysis>                                         \
   using DirichletVertex =                                              \
      ::PACKAGE_NS::detail::DirichletVertex<Analysis, VertexTraits>;	\
   template <typename Analysis>                                         \
   using VertexPointer = std::unique_ptr<Vertex<Analysis> >;

   namespace fi {
   namespace detail {
   template <typename> struct StandardVertexProperties;
   template <typename> struct NeumannVertexProperties;
}

   template <typename> struct VertexTraits;
}

   namespace detail {
   template<typename, template <typename> class> struct Vertex;
   template<typename, template <typename> class> struct StandardVertex;
   template<typename, template <typename> class> struct NeumannVertex;
   template<typename, template <typename> class> struct DirichletVertex;
}

   namespace fi {
   DEFINE_VERTEX_TYPES
}

   //
   // --- Edge types -------------------------------------------------
   //
   namespace fi {
   template <typename> struct GridEdge;
   template <typename> struct StraightGridEdge;
   template <typename> struct CircularGridEdge;

   template <typename Analysis>
   using GridEdgePointer = std::unique_ptr<GridEdge<Analysis> >;

   template <typename Analysis>
   using GridEdgeType = GridEdgePointer<Analysis>;
}

   //
   // --- Graph types ------------------------------------------------
   //
   template <typename> struct GraphProperties;

   namespace detail {
   template <typename,
             template <typename> class,
             template <typename> class > struct GridGraph;
}

   namespace fi {
   //! GridGraph type
   template <typename Analysis>
   using GridGraph = ::PACKAGE_NS::detail::GridGraph<
      Analysis, VertexPointer, GridEdgePointer>;

   //! Pointer to a grid graph
   /*!
    * \tparam Analysis is the analysis type
    *
    * \relates GridGraph
    */
   template<typename Analysis>
   using GridGraphPointer = std::shared_ptr<GridGraph<Analysis> >;
}

   //
   // --- Vertex traits ----------------------------------------------
   //
#define COMMON_VERTEX_TRAITS_TYPES                                      \
   using DataType = typename Analysis::DataType;                        \
   using GridGraphType = GridGraph<Analysis>;                           \
   using VertexDescriptor = std::size_t;                                \
   using EdgeDescriptor =                                               \
      ::boost::detail::csr_edge_descriptor<                             \
                                                VertexDescriptor, std::size_t>; \
   using GridGraphPtr = GridGraphPointer<Analysis>;                     \
   using NeumannVertexProperties =                                      \
      detail::NeumannVertexProperties<Analysis>;                        \
   using StandardVertexProperties =                                     \
      detail::StandardVertexProperties<Analysis>;

   namespace fi {
      template <typename Analysis>
      struct VertexTraits
      {
         COMMON_VERTEX_TRAITS_TYPES
         using EdgeProps = GridEdgePointer<Analysis>;

         using CoefficientQuantity = boost::units::quantity<
            boost::units::si::magnetic_reluctivity, DataType>;
      };
   }
}

#ifdef COMMON_VERTEX_TRAITS_TYPES
#undef COMMON_VERTEX_TRAITS_TYPES
#endif
#ifdef DEFINE_VERTEX_TYPES
#undef DEFINE_VERTEX_TYPES
#endif
#ifdef EXTERN_TEMPLATES
#undef EXTERN_TEMPLATES
#endif

#endif // GRID_GRAPHFWD_HPP_mRjov6k1_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
