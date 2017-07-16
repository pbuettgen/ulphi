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
 *  \brief Linear vertex distribution
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef LINEAR_VERTEX_DISTRIBUTION_HPP_fb8IUiL9_
#define LINEAR_VERTEX_DISTRIBUTION_HPP_fb8IUiL9_

#include "config.h"

#include "utility.hpp"
#include "vertex_distribution.hpp"

namespace PACKAGE_NS {
   //! Linear vertex distribution
   struct LinearVertexDistribution:
      VertexDistribution,
      detail::MakeShared<LinearVertexDistribution>
   {
      //! Initialization
      /*!
       * @param number_of_vertices is the number of vertices
       */
      LinearVertexDistribution(VertexIndex number_of_vertices) :
         VertexDistribution(number_of_vertices)
      {
      }

      //! Calculate the nth vertex' position
      /*!
       * \note The total length is normalized to be one!
       *
       * @param vertex_number is the vertex number for which the
       * position shall be calculated.
       *
       * @return the nth vertex's position
       *
       * @throw std::out_of_range in the case of an invalid vertex number.
       */
      virtual double PositionNthVertex(
         VertexIndex vertex_index) const override;
   };
}

#endif // LINEAR_VERTEX_DISTRIBUTION_HPP_fb8IUiL9_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
