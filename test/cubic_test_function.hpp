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

#ifndef CUBIC_TEST_FUNCTION_HPP_riavFiks_
#define CUBIC_TEST_FUNCTION_HPP_riavFiks_

struct CubicTestFunction
{
   explicit CubicTestFunction(
      double a = 10., double b=-.25, double c=.125, double d=-1./16)
   {
      coeffs_[0]=a;
      coeffs_[1]=b;
      coeffs_[2]=c;
      coeffs_[3]=d;
   }

   template <unsigned kDerivation>
   double Eval(double x);

   double operator() (double x) {
      return Eval<0>(x);
   }

private:
   double coeffs_[4];
};


template <unsigned kDerivation>
double CubicTestFunction::Eval(double x)
{
   static_assert(
      kDerivation<4, "Only derivations up to degree 3 are supported");

   const double kDerivFactors[][4] = {
      { 1., 1., 1., 1. },
      { 0., 1., 2., 3. },
      { 0., 0., 2., 6. },
      { 0., 0., 0., 6. }
   };

   double result = 0.;

   for (std::size_t k=3; kDerivation<k; --k)
   {
      result = (coeffs_[k]*kDerivFactors[kDerivation][k]
		 + result)*x;
   }

   return
      result+coeffs_[kDerivation]*kDerivFactors[kDerivation][kDerivation];
}

#endif
