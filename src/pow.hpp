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
 *  \brief Rise something arithmetical to power
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef POW_HPP_I8Z8pnLx_
#define POW_HPP_I8Z8pnLx_

#include "config.h"

namespace PACKAGE_NS {
   namespace detail {
      //! Helper type for power computation
      /*!
       * Required because only functor templates can be partially specialized
       *
       * \tparam exponent is the exponent
       * \tparam Real is some arithmetic type
       */
      template<unsigned exponent, typename Real>
      struct PowHelper {
         //! Recursively compute the power
         inline auto operator ()(
            Real base) -> decltype(base * PowHelper<exponent - 1, Real>()(base))
         {
            return base * PowHelper<exponent - 1, Real>()(base);
         }
      };

      //! Terminate recursive power computation
      /*!
       * \tparam Real is some arithmetic type
       */
      template<typename Real>
      struct PowHelper<1, Real> {
         Real operator()(Real base)
         {
            return base;
         }
      };
   }

   //! Raise a number to power.
   /*!
    * @tparam exponent is the desired exponent
    * @tparam Real is the floating point type to use
    *
    * @param base is the number to rise to power
    *
    * @return the power
    */
   template<unsigned exponent, typename Real = double>
   inline auto Pow(
      Real base) -> decltype(detail::PowHelper<exponent, Real>()(base))
   {
      return detail::PowHelper<exponent, Real>()(base);
   }
}

#endif // POW_HPP_I8Z8pnLx_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
