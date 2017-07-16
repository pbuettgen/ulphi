#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2012-2016 Philipp Büttgenbach
#
# This file is part of ulphi, a CAE tool and library for computing
# electromagnetic fields.
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
from scipy.fftpack import fft, fftshift
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv
from optparse import OptionParser
from math import pi, sin


### Parameters

# Geometry
domain_depth = Qty(3., "cm")
air_surround = Qty(2., "cm")
coil_dim_x = Qty(10., "mm")
coil_dim_y = Qty(4, "cm")
coil_area = coil_dim_x*coil_dim_y
iron_width = Qty(2., "cm")
coil_space_top = Qty(5., "mm")
coil_space_left_side = Qty(3., "mm")
coil_space_right_side = Qty(7., "mm")
domain_height = coil_dim_y+coil_space_top+iron_width+air_surround
domain_width = (2*iron_width+coil_dim_x+coil_space_left_side
                +coil_space_right_side+air_surround)
grid_division = Qty(12, "1/cm")

# Electro - magnetic
current_density = Qty(2.5, 'MA/m**2')
frequency = Qty("50 Hz")

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


def TaskSetup(time_step):
    """Setup the task."""

    ndx = LinearVertexDistribution(int(domain_width*grid_division)+1)
    ndy = LinearVertexDistribution(int(domain_height*grid_division)+1)

    p0 = (Qty('0m'), Qty('0m'))
    p1 = (domain_width, domain_height)

    grid = GridPatchCartesian(p0, p1, ndx, ndy)

    coil_llc = (iron_width+coil_space_left_side, Qty("0 m"))
    grid.SetLoad(
        c_geo.Rectangle(coil_llc, coil_dim_x, coil_dim_y),
        lambda t: current_density*sin(2*pi*(frequency*t)))

    iron_height_in = coil_dim_y + coil_space_top
    iron_width_total = (2*iron_width+ coil_space_left_side
                        +coil_space_right_side+coil_dim_x)
    iron_rects =[
        c_geo.Rectangle(
            (Qty("0m"), Qty("0m")), iron_width, iron_height_in),
        c_geo.Rectangle(
            (iron_width+coil_dim_x+coil_space_left_side
             +coil_space_right_side, Qty("0m")), iron_width, iron_height_in),
        c_geo.Rectangle(
            (Qty("0m"), iron_height_in), iron_width_total, iron_width)]

    for rect in iron_rects:
        grid.SetMaterial(rect, mat)

    dbc = DirichletCondition()

    grid.west_boundary.SetBoundaryCondition(dbc)
    grid.north_boundary.SetBoundaryCondition(dbc)
    grid.east_boundary.SetBoundaryCondition(dbc)

    grid_compiler = GridCompiler(grid)

    return grid_compiler(time_step)


def VerifyResult():
    """Use a simple model to verify the result."""

    iron_height_in = coil_dim_y + coil_space_top

    average_iron_path_length = (
        2*iron_height_in+2*iron_width+coil_dim_x+coil_space_left_side
        +coil_space_right_side)

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

    points = sp.array( [[x.first.to("mm").getValue(),
                        x.second.to("mm").getValue()]
                        for x in grid_graph.positions])
    potentials = sp.array([x.to("mWb/m").getValue()
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
        [(iron_width+coil_space_left_side).to("mm").magnitude, 0.],
        coil_dim_x.to("mm").magnitude,
        coil_dim_y.to("mm").magnitude, color='red', alpha=.5)
    rect2 = mpatches.Rectangle(
        [0., 0.],
        iron_width.to("mm").magnitude,
        (coil_dim_y+coil_space_top).to("mm").magnitude,
        color='black', alpha=.1)
    rect3 = mpatches.Rectangle(
        [(iron_width+coil_dim_x+coil_space_left_side
          +coil_space_right_side).to("mm").magnitude, .0],
        iron_width.to("mm").magnitude,
        (coil_dim_y + coil_space_top).to("mm").magnitude,
        color='black', alpha=.1)
    rect4 = mpatches.Rectangle(
        [.0, (coil_dim_y+coil_space_top).to("mm").magnitude],
        (2*iron_width+coil_dim_x+coil_space_left_side
         +coil_space_right_side).to("mm").magnitude,
        iron_width.to("mm").magnitude,
        color='black', alpha=.1)

    for r in (rect1, rect2, rect3, rect4):
        ax.add_patch(r)

    plt.axis('equal')
    plt.show()


def PlotPotentialDistribution(path):
    plt.plot(
        [x.first.to("mm").magnitude for x in path.positions],
        [x.to("mWb/m").magnitude for x in path.potentials]
    )
    plt.grid()
    plt.xlabel("x/mm")
    plt.ylabel("A*m/mWb")
    plt.show()


def PlotFluxDensOverTime(time, flux_densities):
    flux_dens = [x.to("T").magnitude for x in flux_densities]
    flux_fft = abs(fft(flux_dens))/len(flux_dens)
    fundamental = flux_fft[1]
    plt.subplot(211)
    x_values = [x.to("ms").magnitude for x in time]
    plt.plot(x_values, flux_dens)
    plt.plot(x_values, [2.*fundamental*sin(2*pi*(frequency*t)) for t in time])
    plt.xlabel("time/ms")
    plt.ylabel("B/T")
    plt.grid()
    plt.subplot(212)
    plt.bar(range(len(flux_fft)), fftshift(flux_fft))
    plt.ylabel("B/T")
    plt.grid()
    plt.show()


def SaveData(time, flux_densities, filename):
    num_values = len(flux_densities)
    flux_dens = sp.array([x.to("T").magnitude for x in flux_densities])
    t_values = sp.array([x.to("ms").magnitude for x in time])
    flux_fft = abs(fft(flux_dens))/num_values
    fundamental = 2*flux_fft[1]*sp.sin(
        2*pi*frequency.to("Hz").magnitude*t_values/1e3)

    csv_writer = csv.writer(
        open(filename + "_time.csv", "w"), dialect="excel-tab")

    csv_writer.writerow(["t", "y1", "y2"])

    for row in sp.c_[t_values, flux_dens, fundamental]:
        csv_writer.writerow(row)

    csv_writer = csv.writer(
        open(filename + "_freq.csv", "w"), dialect="excel-tab")

    for row in sp.c_[range(int(-num_values/2), int(num_values/2)),
                     fftshift(flux_fft)]:
        csv_writer.writerow(row)


if "__main__" == __name__:
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="write report to FILE", metavar="FILE")
    parser.add_option("-p", "--plot",
                      action="store_true", default=False,
                      help="plot the result")

    (options, args) = parser.parse_args()

    number_of_timesteps = 4
    time_step = 1./number_of_timesteps/frequency
    grid_graph = TaskSetup(time_step)

    # Check iron core flux
    path = PathInGridGraph(
        grid_graph, (
            (iron_width+coil_dim_x+coil_space_left_side+coil_space_right_side,
             Qty("0m")),
            (2*iron_width+coil_dim_x+coil_space_left_side+coil_space_right_side,
             Qty("0m"))))

    core_flux = []

    for t in range(number_of_timesteps):
        Newton(grid_graph, 1e-5, 24)
        core_flux.append(path.average_flux_density)
        NextTimeStep(grid_graph)

    time_points = [ x*Qty("1 s") for x in sp.linspace(
        Qty("0 s"), 1./frequency-time_step, number_of_timesteps)]

    if options.plot:
        PlotFluxDensOverTime(time_points, core_flux)

    if options.filename:
        SaveData(time_points, core_flux, options.filename)
