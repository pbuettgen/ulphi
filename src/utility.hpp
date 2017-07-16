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
 *  \brief Useful utilities
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef UTILITY_HPP_NW7X3IJu_
#define UTILITY_HPP_NW7X3IJu_

#include <iterator>
#include <memory>
#include <type_traits>

#include "config.h"

#include <boost/mpl/if.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/utility/value_init.hpp>

namespace PACKAGE_NS {

   //! Determine whether two values are almost equal
   /*!
    *  Determine whether two floating point values are equal within a
    *  given tolerance.
    *
    *  \tparam Float is some floating point type
    *
    *  \tparam kMaxULPs defines the tolerance
    *
    *  \returns true if both values are equal within the given
    *  tolerance.
    */
   template <typename Float, std::size_t kMaxULPs = 4>
   inline bool FloatEqual(Float f1, Float f2) {
      return std::abs(boost::math::float_distance(f1, f2)) <= kMaxULPs;
   }

   //! Determine a value's sign
   template <typename Float>
   inline float Signum(Float x)
   {
      const boost::value_initialized<Float> zero;

      if (get(zero) == x) return .0;

      return x < get(zero) ? -1. : 1.;
   }

   namespace detail {

      //! Check that an iterator is an input iterator
      /*!
       * \tparam Iter is some iterator type
       */
      template <typename Iter>
      using IsInputIterator = std::is_same<
         typename std::iterator_traits<Iter>::iterator_category,
         std::input_iterator_tag>;

      //! Check that an iterator is at least a forward iterator
      /*!
       * \tparam Iter is some iterator type
       */
      template <typename Iter>
      using IsAtLeastForwardIterator = std::is_convertible<
         typename std::iterator_traits<Iter>::iterator_category*,
         std::forward_iterator_tag*>;

      //! Create a new object hold by a shared_ptr
      /*!
       * \tparam T is the type to construct
       */
      template <typename T>
      struct MakeShared
      {
         template<typename ... Args>
         static std::shared_ptr<T> New(Args ... arguments)
         {
            return std::make_shared<T>(arguments...);
         }
      };

      //! Create a new object hold by a unique_ptr
      /*!
       * \tparam T is the type to construct
       */
      template <typename T>
      struct MakeUnique
      {
         template <typename ... Args>
         static std::unique_ptr<T> New(Args ... arguments)
         {
            return std::unique_ptr<T>(new T(arguments ...));
         }
      };

      //! Combine std::add_const and std::add_lvalue_reference
      template<typename T>
      using AddConstAndRef =
         typename std::add_lvalue_reference<
         typename std::add_const<T>::type>::type;

      namespace aux {
         //! Select between pass by value or pass by reference
         /*!
          * If a reference is smaller than the type itself, the reference
          * is prefered.
          *
          * @tparam T is the type about which is to decide how to pass it.
          */
         template<typename T>
         struct PassType {
         private:
            using InType = typename std::remove_reference<T>::type;
         public:
            using OutType = typename boost::mpl::if_c<
            (sizeof(InType)>sizeof(InType&)), const InType&, InType>::type;
         };
      }

      //! Select between pass by value or pass by reference
      template <typename T>
      using PassType = typename aux::PassType<T>::OutType;

      //! Some Property
      /*!
       * \tparam T is the property's type
       *
       * Can be used for object composition.
       */
      template <typename T>
      struct Property
      {
         using TPassType = PassType<T>;

         Property(TPassType t) : t_(t) {}

         void SetProperty(TPassType t)
         {
            t = t_;
         }

         TPassType GetProperty() const
         {
            return t_;
         }

      private:
         T t_;
      };
   }
}

#endif // UTILITY_HPP_NW7X3IJu_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
