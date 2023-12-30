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

from ulphi.materials.models import ArctanReluctivity, SqrtReluctivity, SplineReluctivity
from ulphi.quantities import Qty
from optparse import OptionParser
from sys import exit
import csv
import scipy as sp
import matplotlib.pyplot as plt


m400_50A_H = [ Qty(x, 'A/m') for x in
               [75., 84., 93., 104., 117., 134., 159., 199., 282., 505.,
                1284., 3343., 6787., 11712., 19044]]
m400_50A_J = [ Qty(x, 'T') for x in
               [.5, .6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                 1.7, 1.8, 1.9]]


def PlotReluctivityCharaceristic(h_data, bi_data, model):
    """Plot the material's characteristic."""

    mu0 = Qty(1., "mu0")
    
    brange = [Qty(x, "T") for x in sp.arange(.2, 3.8, .05)]
    nurange = [model(x)[0] for x in brange]

    plt.semilogy(
        [x.inBaseUnits().getValue() for x in brange],
        [y.inBaseUnits().getValue() for y in nurange], label="Interp. Data")

    b_data = [x[1]+mu0*x[0] for x in sp.c_[h_data, bi_data] ]
    nu_data = [ x[0]/(x[1]+mu0*x[0]) for x in sp.c_[h_data, bi_data] ]

    plt.semilogy(
        [x.inBaseUnits().getValue() for x in b_data],
        [y.inBaseUnits().getValue() for y in nu_data],
        "*", label="Orig. Data")

    plt.legend(loc=2)
    plt.ylabel(u"ν*H/m")
    plt.grid()


def WriteOutData(h_data, bi_data, model, outfile):
    """Save data to a file."""
    
    mu0 = Qty(1., "mu0")
    
    b_data = [x[1]+mu0*x[0] for x in sp.c_[h_data, bi_data] ]
    nu_data = [ x[0]/(x[1]+mu0*x[0])for x in sp.c_[h_data, bi_data] ]
    nu_model = [model(x)[0] for x in b_data]
    b_data_plain = [x.inBaseUnits().getValue() for x in b_data]
    nu_data_plain = [x.inBaseUnits().getValue() for x in nu_data]
    nu_model_plain = [x.inBaseUnits().getValue() for x in nu_model]

    csvfile = open(outfile, 'w')
    csvwriter = csv.writer(csvfile, delimiter='\t')
    csvwriter.writerow(["B", "NuOrig", "NuModel"])
    for row in sp.c_[b_data_plain, nu_data_plain, nu_model_plain]:
        csvwriter.writerow(row)
        

if "__main__" == __name__:
    parser = OptionParser()
    parser.add_option(
        "-r", "--reluctivity-model", dest="reluctivity_model",
        help="specify the reluctivity model")
    (options, args) = parser.parse_args()
    
    if "arctan" == options.reluctivity_model:
        reluctivity_model = ArctanReluctivity(m400_50A_H, m400_50A_J)
    elif "sqrt" == options.reluctivity_model:
        reluctivity_model = SqrtReluctivity(m400_50A_H, m400_50A_J)
    elif "spline" == options.reluctivity_model:
        reluctivity_model = SplineReluctivity(m400_50A_H, m400_50A_J)
    else:
        exit("Unknown reluctivity model " + options.reluctivity_model)
        
    WriteOutData(
        m400_50A_H, m400_50A_J, reluctivity_model,
        "M400_50A-"+options.reluctivity_model + ".csv")
