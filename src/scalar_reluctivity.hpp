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
 *  \brief Base class for scalar reluctivity models
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef SCALAR_RELUCTIVITY_HPP_yZ2sPCgc_
#define SCALAR_RELUCTIVITY_HPP_yZ2sPCgc_

#include "config.h"

#include "types.hpp"

namespace PACKAGE_NS {
   //! Base class for scalar magnetic reluctivity models
   struct ScalarReluctivity {

      //! Destruction
      virtual ~ScalarReluctivity();

      //! Calculate the Reluctivity and its derivation
      /*!
       * @param b is the magnetic flux density for which the
       * reluctivity and its derivation shall be computed.
       *
       * @returns a pair of the reluctivity and its derivation.
       */
      virtual ReluctivityReluctivityDerivedPair operator()(
         MagneticFluxDensityReal b) const =0;

      virtual bool ResultCacheRecommended() const =0;
   };

   //! Shared pointer to an immutable reluctivity object
   using ConstReluctivityPointer = std::shared_ptr<const ScalarReluctivity>;
}

#endif // SCALAR_RELUCTIVITY_HPP_yZ2sPCgc_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
