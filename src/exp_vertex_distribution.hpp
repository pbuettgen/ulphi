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
 *  \brief Exponential vertex distribution
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef EXP_VERTEX_DISTRIBUTION_HPP_SNu3dg8H_
#define EXP_VERTEX_DISTRIBUTION_HPP_SNu3dg8H_

#include <cmath>

#include "config.h"

#include <boost/math/constants/constants.hpp>

#include "utility.hpp"
#include "vertex_distribution.hpp"

namespace PACKAGE_NS {
   //! Exponential vertex distribution
   struct ExpVertexDistribution:
      VertexDistribution,
      detail::MakeShared<ExpVertexDistribution>
   {
      //! Initialization
      /*!
       * @param number_of_vertices is the number of vertices
       * @param reverse reverses the distribution
       * @param base is the power's base
       * @throw DomainError in the case of an invalid base.
       */
      ExpVertexDistribution(VertexIndex number_of_vertices, bool reverse,
                            double base);

      //! Compute the n'th vertex's position
      /*!
       * As an exception from the exponential distribution rule the
       * first vertex is located exactly onto the boundary.
       *
       * @param vertex_index is the vertex's index. (Counting starts
       * with zero.)
       *
       * @return the n'th vertex's position
       */
      double PositionNthVertex(VertexIndex vertex_index) const override;

      //! Switch direction
      void Reverse()
      {
         reverse_ = !reverse_;
      }

   private:
      bool reverse_;
      double base_;

      static constexpr double kE = boost::math::double_constants::e;
   };

   //--------------------------------------------------------------
   //
   //   Function implementations
   //
   //--------------------------------------------------------------

   inline ExpVertexDistribution::ExpVertexDistribution(
      VertexIndex number_of_vertices,
      bool reverse = false, double base = kE) :
      VertexDistribution(number_of_vertices), reverse_(reverse), base_(
         std::abs(base))
   {
      if (base_ <= 1.) {
#ifdef __GNUC__
#define FUNCTION_NAME __PRETTY_FUNCTION__
#else
#define FUNCTION_NAME (typeid(*this).name())
#endif
         boost::format error_message(
            "%1%: Invalid base %2%! Base must be greater than one!");
         throw DomainError((error_message % FUNCTION_NAME % base_).str());
#undef FUNCTION_NAME
      }
   }
}

#endif // EXP_VERTEX_DISTRIBUTION_HPP_SNu3dg8H_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */

