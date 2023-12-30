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
from ulphi.geometry import cartesian as c_geo
from ulphi.analysis.magnetostatic import *
from ulphi.materials import Material, ComboReluctivityFluxPar, ConstantReluctivity, SplineReluctivity as NonlinearReluctivity, Axis
import scipy as sp
from scipy import interpolate as ip
from scipy import optimize as opt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv
from optparse import OptionParser


### Parameters

# Geometry
domain_depth = Qty(3., "cm")
coil_dim_x = Qty(5., "mm")
coil_dim_y = Qty(2.4, "cm")
coil_area = coil_dim_x*coil_dim_y
iron_width = Qty(2.5, "cm")
coil_space_top = Qty(5., "mm")
coil_space_side = Qty(5., "mm")
coil1_llc = (coil_space_side, Qty("0 m"))
coil2_llc = (3*coil_space_side + iron_width + coil_dim_x, Qty('0m'))
iron1_llc = (2*coil_space_side+coil_dim_x, Qty('0m'))
iron1_height = coil_dim_y + coil_space_top
iron2_llc = (Qty('0m'), coil_dim_y + coil_space_top)
iron2_width = 2*coil_space_side+coil_dim_x+iron_width
domain_height = coil_dim_y + coil_space_top + 2*iron_width
domain_width = 3*coil_space_side + 2*coil_dim_x + 2*iron_width
grid_division = Qty(20, "1/cm")

# Electro - magnetic
current_density = Qty(1.8, 'MA/m**2')

m400_50A_H = [ Qty(x, 'A/m') for x in
               [75., 84., 93., 104., 117., 134., 159., 199., 282., 505.,
                1284., 3343., 6787., 11712., 19044]]
m400_50A_J = [ Qty(x, 'T') for x in
               [.5, .6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                 1.7, 1.8, 1.9]]
iron_fill_factor = .97
mat = Material(isotropic_reluctivity=ComboReluctivityFluxPar(
        NonlinearReluctivity(m400_50A_H, m400_50A_J),
        ConstantReluctivity(1.), iron_fill_factor))


def PlotMaterialCharaceristic():
    """Plot the material's characteristic."""

    mu0 = Qty(1., "mu0").inBaseUnits()
    nu0 = 1./mu0

    brange = [Qty(x, "T") for x in sp.arange(.2, 3.8, .02)]
    nurange = [mat.EvaluateReluctivity1(x)[0] for x in brange]
    diff_nu = [mat.EvaluateReluctivity1(x)[1] for x in brange]

    plt.subplot(2, 1, 1)
    plt.semilogy(
            [x.getValue() for x in brange],
            [y.getValue() for y in nurange], label="Interp. Data")

    m400_50A_B = [x[1]+mu0*x[0] for x in sp.c_[m400_50A_H, m400_50A_J] ]
    m400_50A_nu = [ x[0]/(x[1]+mu0*x[0])
            for x in sp.c_[m400_50A_H, m400_50A_J] ]
    nu_eff = [ x*nu0/(iron_fill_factor*nu0-iron_fill_factor*x+x)
            for x in m400_50A_nu ]

    plt.semilogy(
            [x.getValue() for x in m400_50A_B],
            [y.getValue() for y in nu_eff], "*", label="Orig. Data")

    plt.legend(loc=2)
    plt.ylabel(u"ν*H/m")
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.semilogy(
            [x.getValue() for x in brange],
            [y.getValue() + 200. for y in diff_nu])
    plt.xlabel("B/T")
    plt.ylabel(u"ν'*H/m + 200")
    plt.grid()
    plt.show()


def TaskSetup():
    """Setup the task."""

    ndx = LinearVertexDistribution(int(domain_width*grid_division)+1)
    ndy = LinearVertexDistribution(int(domain_height*grid_division)+1)

    p0 = (Qty('0m'), Qty('0m'))
    p1 = (domain_width, domain_height)

    grid = GridPatchCartesian(p0, p1, ndx, ndy)

    grid.SetLoad(
        c_geo.Rectangle(coil1_llc, coil_dim_x, coil_dim_y),
        current_density)
    grid.SetLoad(
        c_geo.Rectangle(coil2_llc, coil_dim_x, coil_dim_y),
        -current_density)

    grid.SetMaterial(c_geo.Rectangle(
        iron1_llc, iron_width, iron1_height + .1*iron_width), mat)
    grid.SetMaterial(c_geo.Rectangle(
        iron2_llc, iron2_width, iron_width), mat)

    dbc = DirichletCondition()

    grid.north_boundary.SetBoundaryCondition(dbc)
    grid.east_boundary.SetBoundaryCondition(dbc)

    grid_compiler = GridCompiler(grid)

    return grid_compiler()


def VerifyResult():
    """Use a simple model to verify the result."""

    average_iron_path_length = (
        coil_dim_x+2*coil_space_side+iron_width+coil_dim_y
        +coil_space_top)

    def target_func(b):
        b = Qty(b, "T")
        reluct = mat.EvaluateReluctivity(b, Axis.first)[0]

        return (reluct*b
                -current_density*coil_area/average_iron_path_length
        ).to_base_units().magnitude

    def target_func_deriv(b):
        b = Qty(b, "T")
        rrd = mat.EvaluateReluctivity(b, Axis.first)

        return (rrd[0]+b*rrd[1]).to_base_units().magnitude

    return Qty(opt.newton(target_func, 1., target_func_deriv), "T")


def PlotResult(grid_graph):
    """Plot the result."""

    points = sp.array( [[x[0].to("mm").magnitude,
                        x.second.to("mm").magnitude]
                        for x in grid_graph.positions])
    potentials = sp.array([x.to("mWb/m").magnitude
                           for x in grid_graph.potentials])

    xmin = min(points[:,0])
    xmax = max(points[:,0])
    ymin = min(points[:,1])
    ymax = max(points[:,1])

    grid_x, grid_y = sp.mgrid[
        xmin:xmax:(domain_width*grid_division+1)*1j,
        ymin:ymax:(domain_height*grid_division+1)*1j]

    grid_z = ip.griddata(points, potentials, (grid_x, grid_y), "linear")

    fig, ax = plt.subplots()
    cp = plt.contour(grid_x, grid_y, grid_z)
    plt.clabel(cp, fmt="%g")
    plt.grid()
    plt.xlabel("x/mm")
    plt.ylabel("y/mm")

    rect1 = mpatches.Rectangle(
            [coil1_llc[0].to("mm").magnitude, 0.],
            coil_dim_x.to("mm").magnitude,
            coil_dim_y.to("mm").magnitude, color='red', alpha=.5)
    rect2 = mpatches.Rectangle(
            [coil2_llc[0].to("mm").magnitude, .0],
            coil_dim_x.to("mm").magnitude,
            coil_dim_y.to("mm").magnitude, color='red', alpha=.5)
    rect3 = mpatches.Rectangle(
            [iron1_llc[0].to("mm").magnitude, .0],
            iron_width.to("mm").magnitude,
            (coil_dim_y + coil_space_top).to("mm").magnitude,
            color='black', alpha=.1)
    rect4 = mpatches.Rectangle(
            [x.to("mm").magnitude
                for x in [iron2_llc[0], iron2_llc.second]],
            iron2_width.to("mm").magnitude,
            iron_width.to("mm").magnitude,
            color='black', alpha=.1)

    ap = [ax.add_patch(r) for r in [rect1, rect2, rect3, rect4]]
    plt.axis('equal')
    plt.savefig("SimpleIronCore.svg")
    plt.show()


def PlotPotentialDistribution(path):
    plt.plot(
        [x[0].to("mm").magnitude for x in path.positions],
        [x.to("mWb/m").magnitude for x in path.potentials]
    )
    plt.grid()
    plt.xlabel("x/mm")
    plt.ylabel("A*m/mWb")
    plt.show()


def SaveData(grid_graph, filename):
    csv_writer = csv.writer(open(filename, "w"), dialect="excel-tab")

    csv_writer.writerow(["x", "y", "potential"])

    points = sp.array( [[x[0].to("mm").magnitude,
                         x[1].to("mm").magnitude]
                        for x in grid_graph.positions])
    potentials = sp.array([x.to("mWb/m").magnitude
                           for x in grid_graph.potentials])

    y_old = points[0][1]
    y_new = y_old

    for row in sp.c_[points, potentials]:
        y_old = y_new
        y_new = row[1]
        if y_new != y_old:
            csv_writer.writerow([])
        csv_writer.writerow(row)


if "__main__" == __name__:
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="write report to FILE", metavar="FILE")
    parser.add_option("-p", "--plot",
                      action="store_true", default=False,
                      help="plot the result")

    (options, args) = parser.parse_args()

    grid_graph = TaskSetup()
    Newton(grid_graph, 1e-5, 24)

    # Check iron core flux
    path = PathInGridGraph(
        grid_graph, (iron1_llc, (iron1_llc[0]+iron_width,
                                            iron1_llc[1])))

    flux_dens_simple = VerifyResult()

    print("Average flux density (iron core): " + str(path.average_flux_density))
    print("Average flux density (simple method): " + str(flux_dens_simple))

    if options.plot:
        PlotResult(grid_graph)

    if options.filename:
        SaveData(grid_graph, options.filename)
