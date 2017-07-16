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
import math
import scipy as sp
from scipy.optimize import newton
import matplotlib.pyplot as plt


def SymCalc():
    nu, B, mu0, mur, Js = symbols("nu b mu_0 initial_mu_r_ saturation_polarization_")

    H = nu*B

    nuBexpr = -B + mu0*H + 2*Js/pi*atan(pi*(mur-1)*mu0*H/2/Js)

    diff_nu = - diff(nuBexpr, B)/diff(nuBexpr, nu)

    print(limit(diff_nu, B, 0))


def BH(H):
    Js = 1.7
    mur = 7500
    pi = math.pi
    atan = math.atan
    mu0 = 4e-7*pi

    return mu0*H+2*Js/pi*atan(pi/(2*Js)*(mur-1)*mu0*H)

def nuB(b):
    Js = 1.7
    mur =7500
    mur_1 = mur-1
    pi = math.pi
    atan = math.atan
    mu0 = 4e-7*pi
    atan_arg_fac = pi/(2*Js)*mur_1*mu0

    return newton(
            lambda nu: b*(mu0*nu-1)+2*Js/pi*atan(atan_arg_fac*nu*b),
            50.,
            lambda nu: mu0*b*(1+mur_1/(1+(atan_arg_fac*nu*b)**2)))

H = sp.arange(1e2, 8e5, 100)
B = sp.array([BH(h) for h in H])
nu = H/B

plt.semilogy(B, nu, label="Method 1")

B = sp.arange(.1, 2.8, .1)
nu = sp.array([nuB(b) for b in B])

plt.semilogy(B, nu, "-*", label="Method 2")

plt.grid()
plt.legend(loc=2)
plt.show()
