// -*- coding: utf-8 -*-
/* 
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and library for computing
 * electromagnetic fields.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

/*!
 *  \file
 *
 *  \brief Fill a range with linearly increasing numbers
 *
 *  \author Philipp Büttgenbach
 *
 *  \copyright All rights reserved.
 */


#ifndef LINSPACE_HPP_AmBynUrr_
#define LINSPACE_HPP_AmBynUrr_

#include <iterator>
#include <algorithm>


//! Functor for generating equidistantly spaced values
/*!
 *  \tparam ValueType is some numerical type
 */
template <typename ValueType>
struct LinspaceGenerator
{
   //! Initialize a linspace generator
   /*!
    *  \param xstart marks the start of the desired value range
    *
    *  \param xend marks the end of the desired value range
    *
    *  \param num_values is the desired number of values
    */
   LinspaceGenerator(
      ValueType xstart, ValueType xend, std::size_t num_values)
      : xstart_(xstart), xend_(xend), num_values_(num_values), next_(xstart)
   {}

   //! Compute the next value within the desired range
   /*!
    *  \returns the next value within the desired range
    */
   ValueType operator() () {
      ValueType v = next_;
      next_ += (xend_-xstart_)*(1./static_cast<double>(num_values_-1));
      return v;
   }

private:
   ValueType xstart_, xend_;
   std::size_t num_values_;
   ValueType next_;
};


//! Fill a range with equidistantly spaced values
/*!
 *  \tparam Iterator is some iterator type
 *
 *  \param xbegin marks the start of the desired value range
 *
 *  \param xend marks the end of the desired value range
 *
 *  \param ibegin points to the start of the sequence which shall be
 *  filled with values
 *
 *  \param iend points one step beyond the end of the sequence
 */
template <typename Iterator>
void linspace (
   typename std::iterator_traits<Iterator>::value_type xbegin,
   typename std::iterator_traits<Iterator>::value_type xend,
   Iterator ibegin, Iterator iend)
{
   typedef typename std::iterator_traits<Iterator>::value_type value_type;

   LinspaceGenerator<value_type> gen(
      xbegin, xend, std::distance(ibegin, iend));

   std::generate(ibegin, iend, gen);
}

#endif // LINSPACE_H_
