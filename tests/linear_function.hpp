// -*- coding: utf-8 -*-
/* 
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and
 * library for computing electromagnetic fields.
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
 *  \brief Functors to calculate linear functions
 *
 *  Functors to compute linear functions and their derivations.  The
 *  derivations are implemented by template specialization.
 *
 *  \author Philipp Büttgenbach
 *
 *  \copyright All rights reserved.
 */


#ifndef LINEAR_FUNCTION_HPP_NimpUrsI_
#define LINEAR_FUNCTION_HPP_NimpUrsI_

//! A linear function
/*!
 *  This general template describes all cases with the derivation
 *  higher than degree one.  All other degrees are implemented by
 *  template specialization.
 *
 *  \tparam kDerivation is the derivation's degree.  (Use zero for no
 *  derivation.)
 *
 *  \tparam Floating is some floating point type
 */
template <unsigned kDerivation, typename Floating=double>
struct LinearFunction
{
   //! The functors argument type
   typedef Floating argument_type;

   //! The functors result type
   typedef Floating result_type;

   //! Initialize a linear functor object
   /*!
    *  Does nothing because there is nothing to do in the general
    *  case.
    */
   LinearFunction(Floating, Floating, Floating, Floating)
   {
   }

   //! Calculate the function's value for a specific argument
   /*!
    *  \returns always zero
    */
   result_type operator() (argument_type ) const
   {
      return static_cast<result_type>(0.);
   }
};


//! Specialization for degree zero (no derivation)
/*!
 *  \tparam Floating is some floating point type
 */
template <typename Floating>
struct LinearFunction<0, Floating>
{
   //! The functors argument type
   typedef Floating argument_type;

   //! The functors result type
   typedef Floating result_type;

   //! Initialize a linear functor object
   /*!
    *  This functor object is initialized by giving two points.
    *
    *  \param x0 is the first point's abzissa
    *
    *  \param y0 is the first point's ordinate
    *
    *  \param x1 is the second point's abzissa
    *
    *  \param y1 is the second point's ordinate
    */
   LinearFunction(Floating x0, Floating y0,
		  Floating x1, Floating y1)
      : m_( (y1-y0)/(x1-x0) ), c_(y0-m_*x0)
   {}

   //! Calculate the function's value for a specific argument
   /*!
    *  \param x is the value for which the function shall be evaluated
    *
    *  \returns the functions value
    */
   result_type operator() (argument_type x) const
   {
      return m_*x + c_;
   }

private:
   Floating m_;
   Floating c_;
};


//! Specialization for the first derivation
/*!
 *  \tparam Floating is some floating point type
 */
template <typename Floating>
struct LinearFunction<1, Floating>
{
   //! The functors argument type
   typedef Floating argument_type;

   //! The functors result type
   typedef Floating result_type;

   //! Initialize a linear functor object
   /*!
    *  This functor object is initialized by giving two points.
    *
    *  \param x0 is the first point's abzissa
    *
    *  \param y0 is the first point's ordinate
    *
    *  \param x1 is the second point's abzissa
    *
    *  \param y1 is the second point's ordinate
    */
   LinearFunction(Floating x0, Floating y0,
		  Floating x1, Floating y1)
      : m_( (y1-y0)/(x1-x0) )
   {}

   //! Calculate the function's value for a specific argument
   /*!
    *  The argument is ignored because the first derivation is a
    *  constant.
    *
    *  \returns the first derivation's value
    */
   result_type operator() (argument_type) const
   {
      return m_;
   }

private:
   Floating m_;
};


#endif // LINEARFUNCTION_H_
