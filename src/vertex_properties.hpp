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
 *  \brief Base class for all vertex distributions
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef VERTEX_PROPERTIES_HPP_c9gPc86e_
#define VERTEX_PROPERTIES_HPP_c9gPc86e_

#include "config.h"

#include "load_definition-fwd.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Store a potential
      /*!
       * \tparam Analysis is the analysis type
       *
       * \tparam kNumValues number of values to keep
       */
      template<typename Analysis, std::size_t kNumValues>
      struct PotentialStorage {
      private:
         using MagneticVectorPotential =
            typename Analysis::MagneticVectorPotential;
         using Storage = std::array<MagneticVectorPotential, kNumValues>;

      public:
         //! ConstIterator
         using ConstIterator = typename Storage::const_iterator;

      public:
         //! Initialization
         PotentialStorage()
         {
         }

         //! Initialization
         /*!
          * @param potential initial potential
          */
         explicit PotentialStorage(MagneticVectorPotential potential)
         {
            potentials_[0] = potential;
         }

         //! Get the actual potential
         MagneticVectorPotential Potential() const
         {
            return potentials_[0];
         }

         //! Set potential
         void SetPotential(MagneticVectorPotential new_potential)
         {
            potentials_[0] = new_potential;
         }

         //! Shift stored potential
         /*!
          * The oldest potential value is droped
          */
         void Shift()
         {
            for (std::size_t i = kNumValues - 1; i; --i)
               potentials_[i] = potentials_[i - 1];
         }

         //! Begin
         ConstIterator Begin() const
         {
            return potentials_.begin();
         }

      private:
         Storage potentials_;
      };

      //! Store just one potential value
      template<typename Analysis>
      struct PotentialStorage<Analysis, 1> {
      private:
         using MagneticVectorPotential =
            typename Analysis::MagneticVectorPotential;
         using Storage = MagneticVectorPotential;

      public:
         //! ConstIterator
         using ConstIterator = const MagneticVectorPotential*;

      public:
         //! Default initialization
         PotentialStorage() :
            potential_(.0 * boost::units::si::weber_per_metre)
         {
         }

         //! Initialization
         /*!
          * @param potential initial value
          */
         explicit PotentialStorage(MagneticVectorPotential potential) :
            potential_(potential)
         {
         }

         //! Retrieve potential
         MagneticVectorPotential Potential() const
         {
            return potential_;
         }

         //! Set potential
         void SetPotential(MagneticVectorPotential new_potential)
         {
            potential_ = new_potential;
         }

         //! Shift (no op)
         void Shift()
         {
         }

         //! Begin
         ConstIterator Begin() const
         {
            return &potential_;
         }

      private:
         Storage potential_;
      };

      template<typename Analysis>
      struct PotentialStorage<Analysis, 0> : PotentialStorage<Analysis, 1>
      {
      private:
         using Base = PotentialStorage<Analysis, 1>;
         using MagneticVectorPotential =
            typename Analysis::MagneticVectorPotential;

      public:
         PotentialStorage() :
            Base()
         {
         }

         PotentialStorage(MagneticVectorPotential potential) :
            Base(potential)
         {
         }
      };
   }

   namespace fi {
      namespace detail {
         //! Properties for a standard vertex
         /*!
          *  Used to plug all required properties into a standard vertex.
          *
          *  \tparam Analysis is the analysis type
          */
         template <typename Analysis>
         struct StandardVertexProperties
         {
         private:
            template <typename T>
            using PassType = ::PACKAGE_NS::detail::PassType<T>;

            using DataType = typename Analysis::DataType;
            using CoefficientQuantity = boost::units::quantity<
               boost::units::si::magnetic_reluctivity, DataType>;
            using Current = typename Analysis::Current;
            using LoadDefPointer = ConstLoadDefinitionPointer<Analysis>;
            using MagneticVectorPotential =
               typename Analysis::MagneticVectorPotential;
            using Time = typename Analysis::Time;

         public:
            //! Initialization
            /*!
             *  @param load_definition load definition
             *
             *  @param area secondary grid's cell size
             */
            StandardVertexProperties(
               LoadDefPointer load_definition, Area area)
               : load_definition_(load_definition), area_(area)
            {
            }

         protected:
            //! Contribution to the equation system's right hand side
            Current StandardRhs(Time t) const {
               return area_ * (*load_definition_)(t);
            }

            //! Retrieve potential
            PassType<MagneticVectorPotential> DoGetPotential() const {
               return potential_;
            }

            void DoSetPotential(
               PassType<MagneticVectorPotential> potential)
            {
               potential_ = potential;
            }

         private:
            ConstLoadDefinitionPointer<Analysis> load_definition_;
            Area area_;
            MagneticVectorPotential potential_;
         };
      }
   }

   namespace fi {
      namespace detail {
         //! Properties for a neumann vertex
         /*!
          *  Used to plug all required properties into a neumann vertex.
          *
          *  \tparam Analysis is the analysis type
          */
         template <typename Analysis>
         struct NeumannVertexProperties
         {
         private:
            template <typename T>
            using PassType = ::PACKAGE_NS::detail::PassType<T>;

            using Time = typename Analysis::Time;
            using Current = typename Analysis::Current;
            using MagneticFluxDensity =
               typename Analysis::MagneticFluxDensity;
            using CoefficientQuantity = boost::units::quantity<
               boost::units::si::magnetic_reluctivity,
               typename Analysis::DataType>;

         public:
            //! Default initialization
            NeumannVertexProperties()
               : reluctivity_outside_(.0 * boost::units::si::metre_per_henry),
                 flux_density_outside_(.0 * boost::units::si::tesla),
                 ghost_grid_spacing_(.0 * boost::units::si::metre),
                 parallelity_factor_(Parallelity::kParallel)
            {
            }

            //! Initialization
            /*!
             * @param reluctivity_outside (homogenious) reluctivity
             * outside the computation domain.
             *
             * @param flux_density_outside flux density parallel to
             * the computation domain's boundary.
             *
             * @param delta ghost grid spacing.
             *
             * @param parallelity defines whether the boundary is
             * parallel or antiparallel.
             */
            NeumannVertexProperties(
               MagneticReluctivity reluctivity_outside,
               PassType<MagneticFluxDensity> flux_density_outside,
               Length delta, Parallelity parallelity)
               : reluctivity_outside_(reluctivity_outside),
                 flux_density_outside_(flux_density_outside),
                 ghost_grid_spacing_(delta),
                 parallelity_factor_(parallelity)
            {
            }

         protected:
            //! Contribution to the equation system's right hand side
            Current NeumannRhs() const {
               return static_cast<double>(parallelity_factor_)
                  * reluctivity_outside_
                  * ghost_grid_spacing_
                  * flux_density_outside_;
            }

         private:
            MagneticReluctivity reluctivity_outside_;
            MagneticFluxDensity flux_density_outside_;
            Length ghost_grid_spacing_;
            Parallelity parallelity_factor_;
         };
      }
   }
}

#endif // VERTEX_PROPERTIES_HPP_c9gPc86e_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
