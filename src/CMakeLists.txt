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

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Boost_INCLUDE_DIR}
  ${EIGEN3_INCLUDE_DIR}
  ${LOKI_INCLUDE_DIR})

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library("lib${PROJECT_NAME}_base" OBJECT grid_assembly.cpp)

# Analysis
foreach(analysis_type static harmonic)
  add_library(
    "lib${PROJECT_NAME}_analysis_magneto${analysis_type}" OBJECT
    boundary.cpp
    dirichlet_condition.cpp
    edge_properties.cpp
    grid_compiler.cpp
    grid_graph.cpp
    grid_patch.cpp
    neumann_condition.cpp)
  target_compile_options(
    "lib${PROJECT_NAME}_analysis_magneto${analysis_type}"
    PRIVATE -include magneto${analysis_type}_analysis.hpp)
  target_compile_definitions(
    lib${PROJECT_NAME}_analysis_magneto${analysis_type}
    PRIVATE INSTANTIATE_)
endforeach()

# Geometry
foreach(geometry_type cartesian polar)
  add_library(
    "lib${PROJECT_NAME}_geometry_basic_${geometry_type}"
    OBJECT rectangle.cpp)
  add_library(
    "lib${PROJECT_NAME}_geometry_extra_${geometry_type}"
    OBJECT circle.cpp)
  foreach(geometry_level basic extra)
    target_compile_options(
      lib${PROJECT_NAME}_geometry_${geometry_level}_${geometry_type}
      PRIVATE -include ${geometry_type}_geometry.hpp)
    target_compile_definitions(
      lib${PROJECT_NAME}_geometry_${geometry_level}_${geometry_type}
      PRIVATE INSTANTIATE_)
  endforeach()
endforeach()

# Materials basic
add_library(lib${PROJECT_NAME}_materials_basic OBJECT
  constant_reluctivity.cpp
  scalar_reluctivity.cpp)

# Material models
add_library(lib${PROJECT_NAME}_materials_models OBJECT
  arctan_reluctivity.cpp
  eigen_extra.cpp
  spline_reluctivity.cpp
  sqrt_reluctivity.cpp)

# Materials extra
add_library(lib${PROJECT_NAME}_materials_extra OBJECT
  combined_reluctivity.cpp)

# Vertex distribution
add_library(lib${PROJECT_NAME}_vertex_distributions OBJECT
  exp_vertex_distribution.cpp
  linear_vertex_distribution.cpp
  vertex_distribution.cpp)

# Transformations
add_library(lib${PROJECT_NAME}_transformations OBJECT
  transformation.cpp)

set(part_libs "")
foreach(l
    _analysis_magnetostatic
    _analysis_magnetoharmonic
    _geometry_basic_cartesian
    _geometry_basic_polar
    _geometry_extra_cartesian
    _geometry_extra_polar
    _materials_basic
    _materials_models
    _materials_extra
    _vertex_distributions
    _transformations)
  list(APPEND part_libs lib${PROJECT_NAME}${l})
endforeach()

# Libulphi
set(LIB_ALL_OBJECTS "")
foreach(lib ${part_libs})
  list(APPEND LIB_ALL_OBJECTS $<TARGET_OBJECTS:${lib}>)
endforeach()

add_library(${PROJECT_NAME} SHARED ${LIB_ALL_OBJECTS}
  $<TARGET_OBJECTS:lib${PROJECT_NAME}_base>)
set_target_properties(
  ${PROJECT_NAME} PROPERTIES
  VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}
  SOVERSION ${PROJECT_VERSION_MAJOR})

install(
  TARGETS ${PROJECT_NAME}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  COMPONENT lib)

file(GLOB pkginclude_HEADERS
  LIST_DIRECTORIES false
  ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)

install(
  FILES ${pkginclude_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  COMPONENT dev)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  COMPONENT dev)
