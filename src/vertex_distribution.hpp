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
 *  \brief Base class for all vertex distributions
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef VERTEX_DISTRIBUTION_HPP_YxcpGm7P_
#define VERTEX_DISTRIBUTION_HPP_YxcpGm7P_

#include "config.h"

#include <boost/format.hpp>

#include "exception.hpp"

namespace PACKAGE_NS {

   //! Base class for vertex distributions
   /*!
    * The vertex distribution describes how vertices are distributed
    * along one axis.
    *
    * \note The vertex index starts with zero for the first vertex!
    */
   struct VertexDistribution {
      using VertexIndex = std::size_t;

      //! Initialize a vertex distribution
      /*!
       * @param number_of_vertices is the desired number of vertices within
       * this distribution.
       */
      VertexDistribution(VertexIndex number_of_vertices) :
         number_of_vertices_(number_of_vertices < 2 ? 2 : number_of_vertices)
      {
      }

      //! Destroy a vertex distribution.
      virtual ~VertexDistribution();

      //! Retrieve the number of vertices within this distribution
      /*!
       * @return the number of vertices within this distribution
       */
      std::size_t number_of_vertices() const
      {
         return number_of_vertices_;
      }

      //! Change the number of vertices within this distribution
      /*!
       * Select a new number of vertices within this distribution.
       *
       * @param new_number_of_vertices is the desired number of vertices
       * within this distribution.
       */
      void set_number_of_vertices(VertexIndex new_number_of_vertices)
      {
         number_of_vertices_
            = new_number_of_vertices < 2 ? 2 : new_number_of_vertices;
      }

      //! Calculate the nth vertex's position
      /*!
       * \note The total length is normalized to be one!
       *
       * @param vertex_number is the vertex number for which the position
       * shall be calculated.
       *
       * @return the nth vertex's position
       *
       * @throw std::out_of_range in the case of an invalid vertex number.
       */
      virtual double PositionNthVertex(VertexIndex vertex_index) const =0;

   protected:
      //! Check a vertex number for Validity
      /*!
       * Check whether a vertex number is valid within this vertex distribution.
       *
       * @param vertex_number is the vertex number to check
       *
       * @throw OutOfRange in the case of an invalid vertex number.
       */
      void CheckValid(VertexIndex vertex_index) const
      {
         if (vertex_index >= number_of_vertices_) {
            boost::format error_msg("%1%: Vertex number %2% out of range [0:%3%[!");
#ifdef __GNUC__
#define FUNCTION_NAME __PRETTY_FUNCTION__
#else
#define FUNCTION_NAME (typeid(*this).name())
#endif
            throw OutOfRange(
               (error_msg % FUNCTION_NAME % vertex_index
                % number_of_vertices_).str());
#undef FUNCTION_NAME
         }
      }

   private:
      VertexIndex number_of_vertices_;
   };

   using ConstVertexDistributionPointer =
      std::shared_ptr<const VertexDistribution>;
}

#endif // VERTEX_DISTRIBUTION_HPP_YxcpGm7P_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
