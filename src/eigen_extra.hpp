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
 *  \brief Eigen library addons
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef EIGEN_EXTRA_HPP_CF2XAl1k_
#define EIGEN_EXTRA_HPP_CF2XAl1k_

#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>

#include "pow.hpp"

namespace Eigen {

   //! Compute the arcus tangens from a Eigen matrix or array
   template<typename Array>
   inline auto atan(
      const Array& array) -> decltype(asin(array / sqrt(1. + (array*array))))
   {
      return asin(array / sqrt(1. + array*array));
   }

   //! Check the Levenberg-Marquard-status for an error condition
   /*!
    * @param status is the status which shall be checked whether it
    * indicates an error condition.
    *
    * @throw std::invalid_argument when improper input parameters are
    * detected.
    *
    * @throw std::runtime_error when too many function evaluations are
    * detected or some tolerance is too small.
    */
   void CheckLevenbergMarquardtStatus(
      enum LevenbergMarquardtSpace::Status status);
}

#endif // EIGEN_EXTRA_HPP_CF2XAl1k_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
