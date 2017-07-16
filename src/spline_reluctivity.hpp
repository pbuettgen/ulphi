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
 *  \brief Reluctivity model using a spline
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef SPLINE_RELUCTIVITY_HPP_1FJD0iaY_
#define SPLINE_RELUCTIVITY_HPP_1FJD0iaY_

#include <iterator>
#include <memory>
#include <vector>

#include "config.h"

#include <boost/format.hpp>
#include <boost/utility.hpp>
#include <boost/units/systems/si/codata_constants.hpp>

#include <unsupported/Eigen/Splines>

#include "exception.hpp"
#include "scalar_reluctivity.hpp"
#include "utility.hpp"

namespace PACKAGE_NS {
   //! Model a B-Field dependent magnetic reluctivity using a spline
   struct SplineReluctivity:
      ScalarReluctivity,
      detail::MakeShared<SplineReluctivity>
   {
   private:
      template <typename HIterator, typename BiIterator>
      using AreInputIterators = std::enable_if<
      detail::IsInputIterator<HIterator>::value and
      detail::IsInputIterator<BiIterator>::value>;

      template <typename HIterator, typename BiIterator>
      using AreAtLeastForwardIterators = std::enable_if<
         detail::IsAtLeastForwardIterator<HIterator>::value and
         detail::IsAtLeastForwardIterator<BiIterator>::value>;

   public:
      //! Initialize this reluctivity model from iterators
      /*!
       *  \tparam HIterator is an iterator for specifying a range with
       *  magnetic field intensities
       *
       *  \tparam BiIterator is an iterator for specifying a
       *  range with the corresponding magnetic polarization
       *
       *  \param h_begin should point to the magnetic field intensity
       *  values range's beginning.
       *
       *  \param h_end points just beyond the h-sequence'es end
       *
       *  \param bi_begin should point to the magnetic polarization value
       *  range's beginning.  Make sure that this value range has the
       *  same length as the magnetic field intensity value range!
       *
       */
      template<typename HIterator, typename BiIterator>
      SplineReluctivity(HIterator, HIterator, BiIterator);

      //! Copy construction
      SplineReluctivity(const SplineReluctivity& source)
         : b_min_(source.b_min_), b_max_(source.b_max_),
           spline_(new SplineType(*(source.spline_)))
      {}

   protected:
      //! Do the initialization
      template <typename HIterator, typename BiIterator>
      typename AreInputIterators<HIterator, BiIterator>::type
      Init(HIterator, HIterator, BiIterator);

      //! Do the initialization
      template <typename HIterator, typename BiIterator>
      typename AreAtLeastForwardIterators<HIterator, BiIterator>::type
      Init(HIterator, HIterator, BiIterator);

   public:
      //! Destruct a SplineReluctivity
      ~SplineReluctivity();

      ReluctivityReluctivityDerivedPair operator()(
         MagneticFluxDensityReal b) const override;

      bool ResultCacheRecommended() const override final;

   private:
      using SplineType = Eigen::Spline<double, 1>;
      using SplinePointer = std::unique_ptr<SplineType>;

   private:
      SplineType::Scalar UValue(MagneticFluxDensityReal b) const
      {
         return (b - b_min_) / (b_max_ - b_min_);
      }

   private:
      MagneticFluxDensityReal b_min_, b_max_;
      SplinePointer spline_;
   };

   template<typename HIterator, typename BiIterator>
   SplineReluctivity::SplineReluctivity(
      HIterator h_begin, HIterator h_end, BiIterator bi_begin)
   {
      this->Init(h_begin, h_end, bi_begin);
   }

   template <typename HIterator, typename BiIterator>
   auto SplineReluctivity::Init(
      HIterator h_begin, HIterator h_end, BiIterator bi_begin) ->
      typename AreInputIterators<HIterator, BiIterator>::type
   {
      std::vector<MagneticFieldIntensityReal> h_values(h_begin, h_end);
      std::vector<MagneticPolarizationReal> bi_values;
      const std::size_t number_of_values = h_values.size();
      bi_values.reserve(number_of_values);

      for (std::size_t i=0; i < number_of_values; ++i) {
         bi_values.emplace_back(*bi_begin);
         ++bi_begin;
      }

      this->Init(h_values.begin(), h_values.end(), bi_values.begin());
   }

                                                        template<typename HIterator, typename BiIterator>
                                                        auto SplineReluctivity::Init(
                                                           HIterator h_begin, HIterator h_end, BiIterator bi_begin) ->
                                                           typename AreAtLeastForwardIterators<HIterator, BiIterator>::type
                                                        {
                                                           using boost::units::si::tesla;
                                                           using boost::units::si::ampere_per_metre;
                                                           using boost::units::si::constants::codata::mu_0;

                                                           const std::size_t max_number_of_points = std::distance(h_begin, h_end);

                                                           SplineType::ControlPointVectorType b_values(max_number_of_points),
                                                              nu_values(max_number_of_points);

                                                           SplineType::ControlPointVectorType::Index i = 0;

                                                           while (h_begin != h_end) {
                                                              if ((.0 * tesla != (*bi_begin))
                                                                  and (.0 * ampere_per_metre != (*h_begin))) {
                                                                 b_values(i) = ((*bi_begin) + mu_0 * (*h_begin)) / tesla;
                                                                 nu_values(i) = ((*h_begin) / ampere_per_metre) / b_values(i);
                                                                 ++i;
                                                              }
                                                              ++h_begin;
                                                              ++bi_begin;
                                                           }

                                                           b_min_ = b_values(0) * tesla;
                                                           b_max_ = b_values(i - 1) * tesla;

                                                           if (b_min_ == b_max_) {
                                                              boost::format error_message(
                                                                 "SplineReluctivity: "
                                                                 "Value range seems to be empty: Min: %1%; Max: %2%");
                                                              throw InvalidArgument(
                                                                 (error_message % b_min_.value() % b_max_.value()).str());
                                                           }

                                                           // Normalize B-values
                                                           b_values = (b_values - b_min_ / tesla) / ((b_max_ - b_min_) / tesla);

                                                           spline_ = SplinePointer(
                                                              new SplineType(
                                                                 Eigen::SplineFitting<SplineType>::Interpolate(
                                                                    nu_values.leftCols(i), 3, b_values.leftCols(i))));
                                                        }

                                                                                                                      }

#endif // SPLINE_RELUCTIVITY_HPP_1FJD0iaY_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
