#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2012-2016 Philipp Büttgenbach
#
# This file is part of ulphi, a CAE tool and
# library for computing electromagnetic fields.
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
from ulphi.geometry.polar import *
from ulphi.analysis.magnetostatic import *
from ulphi.materials import Material, ComboReluctivityFluxPar, ComboReluctivityFluxPerp, ConstantReluctivity, SplineReluctivity as NonlinearReluctivity, Axis
from math import pi, sin
import scipy as sp
from scipy import interpolate as ip
from scipy import optimize as opt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv
from optparse import OptionParser


### Parameters

# Geometry
r_in_iron = Qty("9cm")
r_out_iron = Qty("12cm")
r_in_domain = r_in_iron - Qty("5mm")
r_out_domain = r_out_iron + Qty("5mm")
phi_domain = Qty("90deg")
grid_division_r = Qty(24, "1/cm")
grid_division_phi = Qty(.1, "1/deg")

# Electro - magnetic
flux_density = Qty("1.65T")
frequency = Qty("50Hz")
m400_50A_H = [ Qty(x, 'A/m') for x in
               [75., 84., 93., 104., 117., 134., 159., 199., 282., 505.,
                1284., 3343., 6787., 11712., 19044]]
m400_50A_J = [ Qty(x, 'T') for x in
               [.5, .6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                 1.7, 1.8, 1.9]]
m400_50A_reluctivity = NonlinearReluctivity(m400_50A_H, m400_50A_J)
air_reluctivity = ConstantReluctivity(1.)
iron_fill_factor = .97
mat = Material(
        anisotropic_reluctivity=(
            ComboReluctivityFluxPerp(
                m400_50A_reluctivity, air_reluctivity, iron_fill_factor),
            ComboReluctivityFluxPar(
                m400_50A_reluctivity, air_reluctivity, iron_fill_factor)))


def TaskSetup(timestep):
    """Setup the task."""

    vd_r = LinearVertexDistribution(
        int((r_out_domain-r_in_domain)*grid_division_r)+1)
    vd_phi = LinearVertexDistribution(int(phi_domain*grid_division_phi)+1)

    p0 = (r_in_domain, Qty('0deg'))
    p1 = (r_out_domain, phi_domain)

    grid = GridPatchPolar(p0, p1, vd_r, vd_phi)

    grid.SetMaterial(Rectangle(
        (r_in_iron, Qty("0deg")), (r_out_iron, phi_domain)), mat)

    potential_out = flux_density*(r_out_iron-r_in_iron)

    grid.west_boundary.SetBoundaryCondition(DirichletCondition(
        lambda position, time: potential_out*sin(2.*pi*(time*frequency))
    ))
    grid.east_boundary.SetBoundaryCondition(DirichletCondition())

    grid_compiler = GridCompiler(grid)

    return grid_compiler(timestep)


def VerifyResult():
    """Use a simple model to verify the result."""

    average_iron_path_length = (r_in_iron+r_out_iron)/2*phi_domain

    def target_func(b):
        b = Qty(b, "T")
        reluct = mat.EvaluateReluctivity(b, Axis.second)[0]

        return (reluct*b
                -current_density*coil_area/average_iron_path_length
        ).to_base_units().magnitude

    def target_func_deriv(b):
        b = Qty(b, "T")
        rrd = mat.EvaluateReluctivity(b, Axis.second)

        return (rrd[0]+b*rrd[1]).to_base_units().magnitude

    return Qty(opt.newton(target_func, 1., target_func_deriv), "T")


def PlotResult(grid_graph):
    """Plot the result."""

    points = sp.array( [[x[0].to("mm").magnitude,
                        x[1].to("mm").magnitude]
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
                         x.second.to("mm").magnitude]
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

    number_of_timesteps = 16

    grid_graph = TaskSetup(1./frequency/float(number_of_timesteps))

    for k in range(number_of_timesteps):
        Newton(grid_graph, 1e-5, 24)
        NextTimeStep(grid_graph)

    # Check iron core flux
    flux_path = PathInGridGraph(
        grid_graph, ((r_in_domain, Qty("0m")), (r_out_domain, Qty("0m"))))

    r_mid_iron = (r_in_iron + r_out_iron)/2
    mag_volt_path = PathInGridGraph(grid_graph, (
        (r_mid_iron, Qty("0m")), (Qty("0m"), r_mid_iron)))

    if options.plot:
        PlotPotentialDistribution(flux_path)
        #PlotResult(grid_graph)

    if options.filename:
        SaveData(grid_graph, options.filename)
