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
 *  \brief Dirichlet boundary condition
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef DIRICHLET_CONDITION_HPP_3WUgBLnQ_
#define DIRICHLET_CONDITION_HPP_3WUgBLnQ_

#include <functional>

#include "config.h"

#include "analysis.hpp"
#include "boost_units_extra.hpp"
#include "boundary_condition.hpp"
#include "coordinate_systems.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   //! Dirichlet condition
   /*!
    * @tparam Analysis is the analysis type
    */
   template <typename Analysis>
   struct DirichletCondition :
      BoundaryCondition<Analysis>,
      detail::MakeShared<DirichletCondition<Analysis> >
   {
      using ReturnType = typename BoundaryCondition<Analysis>::ReturnType;

      LOKI_DEFINE_CONST_VISITABLE()

      private:
      static constexpr short kBasePriority = 500;

      template <typename T>
      using PassType = detail::PassType<T>;

      using Base = BoundaryCondition<Analysis>;
      using Point = CartesianSystem::Point;
      using Time = typename Analysis::Time;
      using MagneticVectorPotential =
         typename Analysis::MagneticVectorPotential;
      using PotentialFunction =
         std::function<MagneticVectorPotential (const Point&, Time)>;

   public:
      //! Initialization
      /*!
       *  Homogenious Dirichlet condition
       */
      DirichletCondition()
         : Base(kBasePriority), potential_function_(
            [](const Point&, Time) -> MagneticVectorPotential {
               return MagneticVectorPotential::from_value(.0);
            })
      {
      }

      //! Initialization
      /*!
       *  \param fixed_potential a fixed potential
       */
      explicit DirichletCondition(
         PassType<MagneticVectorPotential> fixed_potential)
         : Base(kBasePriority), potential_function_(
            [fixed_potential](const Point&, Time) -> MagneticVectorPotential {
               return fixed_potential;
            })
      {
      }

      //! Set another value
      void set_potential(PassType<MagneticVectorPotential> fixed_potential)
      {
         potential_function_ =
            [fixed_potential](const Point&, Time) -> MagneticVectorPotential {
            return fixed_potential;
         };
      }

      //! Initialization
      /*!
       * @param potential_function takes a point (cartesian) and a
       * time/frequency parameter a returns the related potential.
       */
      explicit DirichletCondition(
         const PotentialFunction& potential_function)
         : Base(kBasePriority), potential_function_(potential_function)
      {
      }

      //! Set a new function for calculating the potential
      void set_potential(const PotentialFunction& potential_function)
      {
         potential_function_ = potential_function;
      }

      //! Retrieve the function for calculating the potential
      const PotentialFunction& potential_function() const {
         return potential_function_;
      }

   private:
      PotentialFunction potential_function_;
   };

#ifndef INSTANTIATE_
   extern template struct DirichletCondition<TimeTransientAnalysis<> >;
   extern template struct DirichletCondition<MagnetoHarmonicAnalysis<> >;
#endif

}

#endif // DIRICHLET_CONDITION_HPP_3WUgBLnQ_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
