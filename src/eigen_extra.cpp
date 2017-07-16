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

#include <stdexcept>

#include "eigen_extra.hpp"

namespace Eigen {
   void CheckLevenbergMarquardtStatus(
      enum LevenbergMarquardtSpace::Status status)
   {
      switch (status) {
         case LevenbergMarquardtSpace::ImproperInputParameters:
            throw std::invalid_argument(
               "Levenberg-Marquardt: Improper input parameters!");
         case LevenbergMarquardtSpace::TooManyFunctionEvaluation:
            throw std::runtime_error(
               "Levenberg-Marquardt: Too many function evaluations!");
         case LevenbergMarquardtSpace::FtolTooSmall:
            throw std::runtime_error(
               "Levenberg-Marquardt: F-tolerance too small!");
         case LevenbergMarquardtSpace::XtolTooSmall:
            throw std::runtime_error(
               "Levenberg-Marquardt: X-tolerance too small!");
         case LevenbergMarquardtSpace::GtolTooSmall:
            throw std::runtime_error(
               "Levenberg-Marquardt: G-tolerance too small!");
         case LevenbergMarquardtSpace::NotStarted:
         case LevenbergMarquardtSpace::Running:
         case LevenbergMarquardtSpace::RelativeReductionTooSmall:
         case LevenbergMarquardtSpace::RelativeErrorTooSmall:
         case LevenbergMarquardtSpace::RelativeErrorAndReductionTooSmall:
         case LevenbergMarquardtSpace::CosinusTooSmall:
         case LevenbergMarquardtSpace::UserAsked:
            ;
      }
   }
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
