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
from ulphi.geometry import polar as p_geo
from ulphi.analysis.magnetostatic import *
from ulphi.materials import *
from math import pi, sin, cos

import scipy as sp
from scipy import interpolate as ip
from scipy import optimize as opt
from scipy.fftpack import fft, fftshift
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from optparse import OptionParser
import csv

# Geometry
iron_r_in = Qty(1., "cm")
iron_width = Qty(2.2, "cm")
iron_len_1 = Qty(5., "cm")
iron_len_2 = Qty(3., "cm")
space_out = 1*iron_width
grid_division = Qty(15, "1/cm")

gp1_llc = (iron_r_in/2., -iron_len_1)
gp1_urc = (iron_r_in+iron_width+space_out, Qty("0m"))
gp2_llc = (iron_r_in/2., Qty("0 rad"))
gp2_urc = (iron_r_in+iron_width+space_out, Qty("90deg"))
gp3_llc = (-iron_len_2, iron_r_in/2)
gp3_urc = (Qty("0m"), iron_r_in+iron_width+space_out)

# Electro - magnetic
peak_flux_density = Qty(1.85, "T")
period = Qty(20., "ms")
frequency = period**(-1)
steps_per_period = 32
time_step = period/steps_per_period
m400_50A_H = [ Qty(x, 'A/m') for x in
               [75., 84., 93., 104., 117., 134., 159., 199., 282., 505.,
                1284., 3343., 6787., 11712., 19044]]
m400_50A_J = [ Qty(x, 'T') for x in
               [.5, .6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                 1.7, 1.8, 1.9]]
iron_fill_factor = .97
m400_50A_reluct = SplineReluctivity(m400_50A_H, m400_50A_J)
air_reluct = ConstantReluctivity(1.)
reluct_flux_par = ComboReluctivityFluxPar(
    m400_50A_reluct, air_reluct, iron_fill_factor)
reluct_flux_perp = ComboReluctivityFluxPerp(
    m400_50A_reluct, air_reluct, iron_fill_factor)
mat1 = Material(anisotropic_reluctivity=(
    reluct_flux_perp, reluct_flux_par))
mat2 = Material(anisotropic_reluctivity=(
    reluct_flux_par, reluct_flux_perp))


def TaskSetup():
    # Grid definitions

    gp1_width = gp1_urc[0] - gp1_llc[0]
    gp1_height = gp1_urc[1] - gp1_llc[1]
    gp1_ndx = LinearVertexDistribution(int(gp1_width*grid_division)+1)
    gp1_ndy = LinearVertexDistribution(int(gp1_height*grid_division)+1)
    gp1 = GridPatchCartesian(gp1_llc, gp1_urc, gp1_ndx, gp1_ndy)

    gp2_nd_phi = LinearVertexDistribution(30)
    gp2 = GridPatchPolar(gp2_llc, gp2_urc, gp1_ndx, gp2_nd_phi)

    gp3_ndx = LinearVertexDistribution(int(iron_len_2*grid_division)+1)
    gp3 = GridPatchCartesian(gp3_llc, gp3_urc, gp3_ndx, gp1_ndx)

    # Grid couplings
    gp2.south_boundary.SetBoundaryCoupling(gp1.north_boundary)
    gp3.east_boundary.SetBoundaryCoupling(gp2.north_boundary)

    # Boundary condition
    bc_dc0 = DirichletCondition()
    gp1.east_boundary.SetBoundaryCondition(bc_dc0)
    gp2.east_boundary.SetBoundaryCondition(bc_dc0)
    gp3.north_boundary.SetBoundaryCondition(bc_dc0)
    gp3.south_boundary.SetBoundaryCondition(bc_dc0)

    # Load definition
    a_amplitude = peak_flux_density * iron_width
    bc_source_field = DirichletCondition(
        lambda pos, t: a_amplitude*cos(2*pi*(frequency*t)))
    gp1.west_boundary.SetBoundaryCondition(bc_source_field)
    gp2.west_boundary.SetBoundaryCondition(bc_source_field)
    gp3.south_boundary.SetBoundaryCondition(bc_source_field)

    # Material definition
    gp1.SetMaterial(
        c_geo.Rectangle(
            (iron_r_in, -iron_len_1), iron_width, iron_len_1), mat1)
    gp2.SetMaterial(
        p_geo.Rectangle(
            (iron_r_in, Qty("0rad")), iron_width, Qty(pi/2, "rad")), mat1)
    gp3.SetMaterial(
        c_geo.Rectangle(
            (-iron_len_2, iron_r_in), iron_len_2*1.1, iron_width), mat2)

    grid_assembly = GridAssembly([gp1, gp2, gp3])
    grid_compiler = GridCompiler(grid_assembly)

    return grid_compiler(time_step)


def PlotPotentialDistribution(path):
    plt.plot(
        [x[0].to("mm").magnitude for x in path.positions],
        [x.to("mWb/m").magnitude for x in path.potentials]
    )
    plt.grid()
    plt.xlabel("x/mm")
    plt.ylabel("A*m/mWb")
    plt.show()


def PlotResult(grid_graph):
    len_unit = "cm"
    pot_unit = "mWb/m"

    points = sp.array([[x[0].to(len_unit).magnitude,
                        x[1].to(len_unit).magnitude]
                       for x in grid_graph.positions])
    potentials = sp.array([x.to(pot_unit).magnitude
                           for x in grid_graph.potentials])

    xmin = min(points[:,0])
    xmax = max(points[:,0])
    ymin = min(points[:,1])
    ymax = max(points[:,1])

    domain_width = iron_len_2+iron_r_in+iron_width+space_out
    domain_height = iron_len_1+iron_r_in+iron_width+space_out

    grid_x, grid_y = sp.mgrid[
        xmin:xmax:domain_width*grid_division*1j,
        ymin:ymax:domain_height*grid_division*1j]

    grid_z = ip.griddata(points, potentials, (grid_x, grid_y), "linear")

    fig, ax = plt.subplots()
    cp = plt.contour(grid_x, grid_y, grid_z)
    plt.clabel(cp)
    plt.grid()
    plt.axis("equal")
    plt.show()


def PlotMagVolt(time_points, mag_voltage, mag_volt_base):
    """Plot magnetic voltage over time"""

    time_unit = "ms"
    mag_volt_unit = "kA"

    x = [ t.to(time_unit).magnitude for t in time_points ]
    y1 = [ v.to(mag_volt_unit).magnitude for v in mag_voltage ]
    y2 = [ v.to(mag_volt_unit).magnitude for v in mag_volt_base ]

    plt.plot(x, y1, "o-")
    plt.plot(x, y2)
    plt.xlabel("Time t/"+time_unit)
    plt.ylabel("Mag. voltage Um/"+mag_volt_unit)
    plt.grid()
    plt.show()


def SaveData(filename, time, mag_volt, mag_volt_base, mv_harm_order, mv_fft):
    time = [x.to("ms").magnitude for x in time]
    mag_volt = [x.to("kA").magnitude for x in mag_volt]
    mag_volt_base = [
        x.to("kA").magnitude for x in mag_volt_base]
    mag_volt_fft = [ x.to("kA").magnitude for x in mv_fft]

    csv_writer = csv.writer(
        open(filename + "_time.csv", "w"), dialect="excel-tab")

    csv_writer.writerow(["t", "y1", "y2"])

    for k in range(len(time)):
        csv_writer.writerow([time[k], mag_volt[k], mag_volt_base[k]])

    csv_writer = csv.writer(
        open(filename + "_freq.csv", "w"), dialect="excel-tab")

    for k in range(len(mv_harm_order)):
        csv_writer.writerow([mv_harm_order[k], mag_volt_fft[k]])


if "__main__" == __name__:
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="write report to FILE", metavar="FILE")
    parser.add_option("-p", "--plot",
                      action="store_true", default=False,
                      help="plot the result")

    (options, args) = parser.parse_args()

    grid_graph = TaskSetup()

    mag_voltage_path = PathInGridGraph(
        grid_graph, (
            (iron_r_in+iron_width/2., -iron_len_1),
            (iron_r_in+iron_width/2., Qty("0m")),
            (Qty("0m"), iron_r_in+iron_width/2.),
            (-iron_len_2, iron_r_in+iron_width/2.)))

    mag_volts = []
    time_points = [time_step*k for k in range(steps_per_period)]

    for k in range(int(steps_per_period/4)):
        Newton(grid_graph, 1e-5, 24)
        mag_volts.append(mag_voltage_path.magnetic_voltage)
        NextTimeStep(grid_graph)

    # Construct a whole cycle from the quater cycle
    mag_volts.append(Qty(.0, 'A'))
    mag_volts.extend([x*(-1) for x in mag_volts[-2::-1]])
    mag_volts.extend(mag_volts[-2:0:-1])

    # Fourier analysis
    mag_volts_fft = [
        abs(x)/steps_per_period for x in
        fft([ y.to_base_units().magnitude for y in mag_volts])]

    mag_volt_base = [
        mag_volts_fft[1]*2.*cos(2*pi*(frequency*t)) for t in time_points]

    mag_volts_fft = [Qty(x, "A") for x in fftshift(mag_volts_fft)]
    mag_volt_base = [Qty(x, "A") for x in mag_volt_base]
    mag_volts_fft = [Qty(x, "A") for x in mag_volts_fft]
    mag_volt_harm_order = range(-int(steps_per_period/2),
                                int(steps_per_period/2))

    if options.plot:
        PlotMagVolt(time_points, mag_volts, mag_volt_base)

    # Save result
    if options.filename:
        SaveData(options.filename,
                 time_points, mag_volts, mag_volt_base,
                 mag_volt_harm_order, mag_volts_fft)
