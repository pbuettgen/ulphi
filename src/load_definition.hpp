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
 *  \brief load definition
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef LOAD_DEFINITION_HPP_0T9OlJQ3_
#define LOAD_DEFINITION_HPP_0T9OlJQ3_

#include <functional>

#include "config.h"

#include "load_definition-fwd.hpp"
#include "shape.hpp"
#include "types.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {

   //! Load definition
   /*!
    * \tparam Analysis is the analysis type
    */
   template <typename Analysis>
   struct LoadDefinition : detail::MakeShared<LoadDefinition<Analysis> >
   {
   private:
      using Time = typename Analysis::Time;
      using Current = typename Analysis::Current;
      using CurrentFunction = std::function<Current (Time)>;
      using CurrentDensity = typename Analysis::CurrentDensity;
      using CurrentDensityFunction = std::function<CurrentDensity (Time)>;

   public:
      //! Initialization
      /*!
       * Constant (time invariant) current density
       *
       * @param current_density current density
       */
      LoadDefinition(CurrentDensity current_density)
	 : current_density_function_(
	    [current_density](Time){ return current_density; })
      {
      }

      //! Initialization
      /*!
       * Constant current homogeniously distributed over a given area.
       *
       * @param current current.
       *
       * @param area area.
       */
      LoadDefinition(Current current, Area area)
	 : current_density_function_(
	    [current, area](Time) -> CurrentDensity { return current/area; })
      {
      }

      //! Initialization
      /*!
       * Time/frequency variing current density
       *
       * @param current_density_function takes the time/frequency as
       * parameter and computes the related current density.
       */
      LoadDefinition(
	 CurrentDensityFunction current_density_function)
	 : current_density_function_(current_density_function)
      {
      }

      //! Initialization
      /*!
       * Time/frequency variing current
       *
       * @param current_function takes the time/frequency as parameter
       * and computes the related current.
       */
      LoadDefinition(CurrentFunction current_function, Area area)
	 : current_density_function_(
	    [current_function, area](Time t) -> CurrentDensity {
	       return current_function(t)/area;
	    })
      {
      }

      //! Compute load
      CurrentDensity operator() (Time t) const
      {
	 return current_density_function_(t);
      }

   private:
      CurrentDensityFunction current_density_function_;
   };
}

#endif // LOAD_DEFINITION_HPP_0T9OlJQ3_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
