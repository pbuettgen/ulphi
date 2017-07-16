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
 *  \brief Index type.
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef INDICES_HPP_Gv98F2i0_
#define INDICES_HPP_Gv98F2i0_

#include <iostream>

#include "config.h"

#include "utility.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Column index tag
      struct ColumnIndexTag {};
      //! Row index tag
      struct RowIndexTag {};
      //! Boundary index tag
      struct BoundaryIndexTag {};

#define UNARY_OP_INDEX(OP)                      \
      Index& operator OP (const Index other) {	\
         i_ OP other.i_;                        \
         return *this;                          \
      }

#define UNARY_OP_SIZE_T(OP)                     \
      Index& operator OP (const size_t other) {	\
         i_ OP other;                           \
         return *this;                          \
      }

#define UNARY_OP_COMPARE(OP)                        \
      bool operator OP (const Index other) const {	\
         return i_ OP other.i_;                     \
      }

#define UNARY_OP_COMPARE_SIZE_T(OP)                 \
      bool operator OP (const size_t other) const {	\
         return i_ OP other;                        \
      }

      //! Generic index type
      /*!
       * \tparam Tag is a tag to make this type unique
       */
      template <typename Tag>
      struct Index
      {
         using TagType = Tag;

         //! size type
         using size_t = std::size_t;

         //! Initialization
         explicit Index(size_t i) : i_(i)
         {
         }

         UNARY_OP_INDEX(+=)
         UNARY_OP_INDEX(-=)
         UNARY_OP_INDEX(*=)
         UNARY_OP_INDEX(=)
         UNARY_OP_SIZE_T(+=)
         UNARY_OP_SIZE_T(-=)
         UNARY_OP_SIZE_T(*=)
         UNARY_OP_COMPARE(==)
         UNARY_OP_COMPARE(<=)
         UNARY_OP_COMPARE(<)
         UNARY_OP_COMPARE(>=)
         UNARY_OP_COMPARE(>)
         UNARY_OP_COMPARE_SIZE_T(==)
         UNARY_OP_COMPARE_SIZE_T(<=)
         UNARY_OP_COMPARE_SIZE_T(<)
         UNARY_OP_COMPARE_SIZE_T(>=)
         UNARY_OP_COMPARE_SIZE_T(>)

         //! Increment
         Index& operator++() {
            ++i_;
            return *this;
         }

         //! Decrement
         Index& operator--() {
            --i_;
            return *this;
         }

         //! Access the bare value
         size_t Get() const {
            return i_;
         }

      private:
         size_t i_;
      };

#define DEFINE_OPERATOR(OP, LHS_TYPE, RHS_TYPE)                     \
      template <typename Tag>                                       \
      LHS_TYPE operator OP (const LHS_TYPE lhs, const RHS_TYPE rhs)	\
      {                                                             \
         LHS_TYPE tmp(lhs);                                         \
         tmp OP ## = rhs;                                           \
         return tmp;                                                \
      }

      DEFINE_OPERATOR(+, Index<Tag>, Index<Tag>)
      DEFINE_OPERATOR(+, Index<Tag>, std::size_t)
      DEFINE_OPERATOR(-, Index<Tag>, Index<Tag>)
      DEFINE_OPERATOR(-, Index<Tag>, std::size_t)
      DEFINE_OPERATOR(*, Index<Tag>, Index<Tag>)
      DEFINE_OPERATOR(*, Index<Tag>, std::size_t)

      //! Output to a stream
      template<typename Tag, typename Character>
      inline std::basic_ostream<Character>& operator <<(
         std::basic_ostream<Character>& ostream, const Index<Tag>& index)
      {
         ostream << index.Get();
         return ostream;
      }

      //! Comparison
      template <typename Tag>
      inline bool operator == (std::size_t lhs, Index<Tag> rhs) {
         return rhs == lhs;
      }
   }

   //! Row index type
   using RowIndex = detail::Index<detail::RowIndexTag>;

   //! Column index type
   using ColumnIndex = detail::Index<detail::ColumnIndexTag>;

   //! Boundary index type
   using BoundaryIndex = detail::Index<detail::BoundaryIndexTag>;

#ifdef UNARY_OP_INDEX
#undef UNARY_OP_INDEX
#endif
#ifdef UNARY_OP_SIZE_T
#undef UNARY_OP_SIZE_T
#endif
#ifdef UNARY_OP_COMPARE
#undef UNARY_OP_COMPARE
#endif
#ifdef UNARY_OP_COMPARE_SIZE_T
#undef UNARY_OP_COMPARE_SIZE_T
#endif
#ifdef DEFINE_OPERATOR
#undef DEFINE_OPERATOR
#endif
}

#endif // INDICES_HPP_Gv98F2i0_
