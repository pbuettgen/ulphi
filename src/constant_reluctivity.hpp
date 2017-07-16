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
 *  \brief Constant reluctivity
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef CONSTANT_RELUCTIVITY_HPP_SpQyWJ9X_
#define CONSTANT_RELUCTIVITY_HPP_SpQyWJ9X_

#include "config.h"

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/codata_constants.hpp>

#include "scalar_reluctivity.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
      //! Reluctivity which is independent from the magnetic field
      struct ConstantReluctivity:
            ScalarReluctivity, detail::MakeShared<ConstantReluctivity>
      {

            //! Construct a ConstantReluctivity from the relative permeability
            /*!
             * The material's absolute reluctivity is calculated from the
             * given relative permeability.
             *
             * @param mu_r is the material's relative permeability
             */
            explicit ConstantReluctivity(double mu_r) :
                  nu_(1. / mu_r / boost::units::si::constants::codata::mu_0)
            {
            }

            //! Construct a ConstantReluctivity from the absolute permeability
            /*!
             * The material's absolute reluctivity is calculated from the
             * given permeability.
             *
             * @param mu is the material's absolute permeability.
             */
            explicit ConstantReluctivity(Permeability mu) :
                  nu_(1. / mu)
            {
            }

            //! Take simply the material's reluctivity
            /*!
             * @param nu is the material's reluctivity
             */
            explicit ConstantReluctivity(MagneticReluctivity nu) :
                  nu_(nu)
            {
            }

            //! Compute reluctivity
            /*!
             * @param b is the magnetic flux density for which the
             * reluctivity shall be computed
             *
             * @return reluctivity and its derivation
             */
            ReluctivityReluctivityDerivedPair operator()(
                  MagneticFluxDensityReal b) const override;

            bool ResultCacheRecommended() const override final;

      private:
            MagneticReluctivity nu_;
      };
}


#endif // CONSTANT_RELUCTIVITY_HPP_SpQyWJ9X_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
