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
 *  \brief Neumann boundary condition
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef NEUMANN_CONDITION_HPP_nTOJwEi3_
#define NEUMANN_CONDITION_HPP_nTOJwEi3_

#include "config.h"

#include <boost/units/systems/si/codata/universal_constants.hpp>

#include "analysis.hpp"
#include "boundary.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   //! Neumann boundary condition
   /*!
    * \tparam Analysis is the analysis type
    */
   template<typename Analysis>
   struct NeumannCondition:
      BoundaryCondition<Analysis>,
      detail::MakeShared<NeumannCondition<Analysis> >
   {
   private:
      using MagneticFluxDensity = typename Analysis::MagneticFluxDensity;
      using Base = BoundaryCondition<Analysis>;
      using ReturnType = typename Base::ReturnType;
      static constexpr short kBasePriority = 50;

   public:
      LOKI_DEFINE_CONST_VISITABLE()

      //! Initialization
      /*!
       * @param reluctivity_outside is the reluctivity outside the
       * finite difference domain.
       *
       * @param flux_density_outside is the * flux density outside the
       * finite difference domain.
       */
      NeumannCondition(MagneticReluctivity reluctivity_outside,
                       MagneticFluxDensity flux_density_outside)
      : Base(kBasePriority + static_cast<short>(
                std::abs(in_base_units(flux_density_outside)))),
         reluctivity_outside_(reluctivity_outside),
         flux_density_outside_(flux_density_outside)
      {
      }

      //! Initialization
      /*!
       * @param permeability_outside is the permeability outside the
       * finite difference domain.
       *
       * @param flux_density_outside is the flux density outside the
       * finite difference domain.
       */
      NeumannCondition(Permeability permeability_outside,
                       MagneticFluxDensity flux_density_outside)
         : Base(kBasePriority + static_cast<short>(
                   std::abs(in_base_units(flux_density_outside)))),
           reluctivity_outside_(1. / permeability_outside),
           flux_density_outside_(flux_density_outside)
      {
      }

      NeumannCondition(double rel_permeability_outside,
                       MagneticFluxDensity flux_density_outside)
         : Base(kBasePriority + static_cast<short>(
                   std::abs(in_base_units(flux_density_outside)))),
           reluctivity_outside_(
              1. / (boost::units::si::constants::codata::mu_0
                    * rel_permeability_outside)),
           flux_density_outside_(flux_density_outside)
      {
      }

      MagneticReluctivity reluctivity_outside() const {
         return reluctivity_outside_;
      }

      MagneticFluxDensity flux_density_outside() const {
         return flux_density_outside_;
      }

   private:
      MagneticReluctivity reluctivity_outside_;
      MagneticFluxDensity flux_density_outside_;
   };

#ifndef INSTANTIATE_
   extern template struct NeumannCondition<TimeTransientAnalysis<> >;
   extern template struct NeumannCondition<MagnetoHarmonicAnalysis<> >;
#endif
}

#endif // NEUMANN_CONDITION_HPP_nTOJwEi3_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
