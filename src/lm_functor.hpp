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
 *  \brief Functor for use with the Levenberg-Marquardt-Procedure
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef LM_FUNCTOR_HPP_rA1EPpmD_
#define LM_FUNCTOR_HPP_rA1EPpmD_

#include <iterator>

#include "config.h"

#include <Eigen/Core>

#include <boost/units/systems/si.hpp>
#include <boost/utility.hpp>

#include "boost_units_extra.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Functor required by the Levenberg-Marquardt procedure
      /*!
       * @tparam NUMBER_OF_PARAMETERS is the number of searched
       * parameters in the target function.
       */
      template<std::size_t NUMBER_OF_PARAMETERS>
      struct LMFunctor {
      private:
         template<typename T>
         using AddConstAndRef = detail::AddConstAndRef<T>;

         template <typename HIterator, typename BiIterator>
         using AreInputIterators = std::enable_if<
            detail::IsInputIterator<HIterator>::value and
            detail::IsInputIterator<BiIterator>::value, int>;

         template <typename HIterator, typename BiIterator>
         using AreAtLeastForwardIterators = std::enable_if<
            detail::IsAtLeastForwardIterator<HIterator>::value and
            detail::IsAtLeastForwardIterator<BiIterator>::value, int>;

      protected:
         //! Array for storing the b_i/h-value-pairs
         using ValueArray = Eigen::Array<double, Eigen::Dynamic, 2>;
         //! Index type for this array
         using ArrayIndex = typename ValueArray::Index;

      public:
         //! Initialization
         /*!
          * @tparam HIterator is a type for iterating over magnetic
          * field intensity values
          *
          * @tparam BiIterator is a type for iterating over magnetic
          * polarization values
          *
          * @param h_begin points to the start of magnetic field
          * intensity values.
          *
          * @param h_end points just beyond the end of these values
          *
          * @param bi_begin points to the start of the corresponding
          * polarization values.
          */
         template<typename HIterator, typename BiIterator,
                  typename AreAtLeastForwardIterators<HIterator,
                                                      BiIterator>::type =0>
         LMFunctor(
            HIterator h_begin, AddConstAndRef<HIterator> h_end,
            BiIterator bi_begin)
            : bi_h_values_(std::distance(h_begin, h_end), 2)
         {
            this->CopyValues(h_begin, h_end, bi_begin);
         }

      private:
         template <typename HIterator, typename BiIterator>
         void CopyValues(
            HIterator h_begin, AddConstAndRef<HIterator> h_end,
            BiIterator bi_begin)
         {
            for (ArrayIndex i = 0; i < bi_h_values_.rows();
                 ++i, ++h_begin, ++bi_begin) {
               bi_h_values_(i, 0) = in_base_units(*h_begin);
               bi_h_values_(i, 1) = in_base_units(*bi_begin);
            }
         }

      public:
         template <typename HIterator, typename BiIterator,
                   typename AreInputIterators<HIterator, BiIterator>::type =0 >
         LMFunctor(
            HIterator h_begin, AddConstAndRef<HIterator> h_end,
            BiIterator bi_begin)
         {
            using MagFieldIntensity
               = typename std::iterator_traits<HIterator>::value_type;
            using MagPolarization
               = typename std::iterator_traits<BiIterator>::value_type;
            std::vector<MagFieldIntensity> h_values(h_begin, h_end);
            std::vector<MagPolarization> bi_values;
            std::size_t number_of_values = h_values.size();
            bi_values.reserve(number_of_values);

            for(std::size_t i=0; i < number_of_values; ++i) {
               bi_values.emplace_back(*bi_begin);
               ++bi_begin;
            }

            bi_h_values_.resize(number_of_values, 2);

            this->CopyValues(
               h_values.begin(), h_values.end(), bi_values.begin());
         }

         //! Number of parameters
         static constexpr int inputs()
         {
            return NUMBER_OF_PARAMETERS;
         }

         //! Number of data points
         int values()
         {
            return bi_h_values_.rows();
         }

      protected:
         ValueArray bi_h_values_;
      };
   }
}

#endif // LM_FUNCTOR_HPP_rA1EPpmD_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
