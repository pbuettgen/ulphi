#! /usr/bin/env python3
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

from ulphi import Qty
from ulphi.geometry.polar import Rectangle as PRectangle
from ulphi.analysis.magnetostatic import *
import scipy as sp
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
import matplotlib.pyplot as plt
from math import pi, log
import csv
from optparse import OptionParser


# outer radius inner conductor
r1 = Qty(5., "mm")
# inner radius outer conductor
r2 = 2.5*r1
# outer radius outer conductor
r3 = 3.5*r1
# outer domain radius
r4 = 4*r1
# angle
domain_phi = Qty(pi/8, "rad")

current = Qty(16, "A")


def PlotResult(x, y1, y2):
    plt.plot(x*1e3, y1*1e6, label="num. solution")
    plt.plot(x*1e3, y2*1e6, label="sym. solution")
    plt.xlabel("Radius r/mm")
    plt.ylabel("Vector potential A*Mm/Wb")
    plt.legend()
    plt.grid()
    plt.show()


def SaveData(filename, data, header):
    csv_file = file(filename, "w")
    csv_writer = csv.writer(csv_file, dialect="excel-tab")

    csv_writer.writerow(header)

    for row in data:
        csv_writer.writerow(row)


def AnalyticalSolution(radius):
    I = Qty(2*pi, "rad")/domain_phi*current
    J1 = I/r1/r1/pi
    mu0 = Qty(1., "mu_0")

    A1 = lambda r: mu0*J1/4.*r*r
    A2 = lambda r: mu0*J1*r1*r1/2*(.5+log(r/r1))
    A3 = lambda r: mu0*J1*r1*r1/2*(
        .5+log(r2/r1)+log(r/r2)*(1+r2*r2/(r3*r3-r2*r2))+(
        r2*r2-r*r)/(2.*(r3*r3-r2*r2)))

    A0 = A3(r3)

    if radius <= r1:
        return A0-A1(radius)
    elif radius > r1 and radius <= r2:
        return A0-A2(radius)
    elif radius > r2 and radius <= r3:
        return A0-A3(radius)
    else:
        return A0-A3(r3)

def NumericalSolution():
    nd_r = LinearVertexDistribution(51)
    nd_phi = LinearVertexDistribution(13)

    p1 = (r1/(nd_r.number_of_vertices-1), Qty("0rad"))
    p2 = (r1, domain_phi)
    gp1 = GridPatchPolar(p1, p2, nd_r, nd_phi)
    gp1.SetLoad(PRectangle(p1, p2), current)

    gp2 = GridPatchPolar(
        (r1, Qty("0rad")), (r2, domain_phi), nd_r, nd_phi)

    p3 = (r2, Qty("0rad"))
    p4 = (r3, domain_phi)
    gp3 = GridPatchPolar(p3, p4, nd_r, nd_phi)
    gp3.SetLoad(PRectangle(p3, p4), -current)

    gp4 = GridPatchPolar(
        (r3, Qty("0rad")), (r4, domain_phi), nd_r, nd_phi)
    gp4.east_boundary.SetBoundaryCondition(DirichletCondition())

    gp1.east_boundary.SetBoundaryCoupling(gp2.west_boundary)
    gp2.east_boundary.SetBoundaryCoupling(gp3.west_boundary)
    gp3.east_boundary.SetBoundaryCoupling(gp4.west_boundary)

    grid_assembly = GridAssembly( [gp1, gp2, gp3, gp4] )
    grid_compiler = GridCompiler(grid_assembly)
    grid_graph = grid_compiler()
    LinearSolve(grid_graph)

    return PathInGridGraph(
        grid_graph, ((r1/10., Qty("0m")), (r4, Qty("0m"))))


def main(options, args):
    path = NumericalSolution()

    path_x = sp.array(
        [x[0].to("m").magnitude for x in path.positions])
    num_solution = sp.array(
        [x.to("Wb/m").magnitude for x in path.potentials])
    interpolator = Spline(path_x, num_solution)
    x_values = sp.linspace(path_x[0], path_x[-1], 13)
    y_values = interpolator(x_values)
    sym_solution = sp.array(
        [AnalyticalSolution(Qty(x, "m")).to("Wb/m").magnitude
         for x in x_values])

    if options.plot:
        PlotResult(x_values, y_values, sym_solution)

    if options.filename:
        SaveData(
            options.filename,
            sp.c_[1e3*x_values, 1e6*interpolator(x_values), 1e6*sym_solution],
             ["x", "NumA", "SymA"])


if "__main__" == __name__:
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="write report to FILE", metavar="FILE")
    parser.add_option("-p", "--plot",
                      action="store_true", default=False,
                      help="plot the result")

    (options, args) = parser.parse_args()

    main(options, args)
