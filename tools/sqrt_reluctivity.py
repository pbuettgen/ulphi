#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2012-2016 Philipp Büttgenbach
#
# This file is part of ulphi, a CAE tool and library for
# computing electromagnetic fields.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#

from __future__ import print_function, division
from sympy import *
from sympy.printing import ccode

nu, B, mu0, mur, Js, a = symbols("nu b mu_0 initial_mu_r_ saturation_polarization_ a_ ")

H = nu*B

Ha = mu0*H*(mur-1)/Js

nuBexpr = -B + mu0*H +Js*(Ha+1-sqrt((Ha+1)**2-4*Ha*(1-a)))/(2*(1-a))

nu_explicit = solve(nuBexpr, nu)

for solution in nu_explicit:
    print(solution.subs({mu0:4e-7, B:1., mur:800., Js:1.8, a:.2}))

print("Expresion for nu:")
print(ccode((nu_explicit[0]), assign_to="nu"))

print("Limit nu: B->0")
print(limit(nu_explicit[0], B, 0))

diff_nu = diff(nu_explicit[0], B)
print("Derivation dnu/dB:")
print(ccode((diff_nu), assign_to="nu_diff"))
print("Derivation's limit B->0")
print(limit(diff_nu, B, 0))
