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

foreach(geometry_type cartesian polar)
  PYTHON_ADD_MODULE(
    _${geometry_type} geometry.cpp
    $<TARGET_OBJECTS:lib${PROJECT_NAME}_geometry_basic_${geometry_type}>
    $<TARGET_OBJECTS:lib${PROJECT_NAME}_geometry_extra_${geometry_type}>)
  target_compile_definitions(
    _${geometry_type} PRIVATE MODULE_NAME=_${geometry_type})
  target_compile_options(
    _${geometry_type} PRIVATE -include ${geometry_type}_geometry.hpp)
endforeach()

install(
  TARGETS _cartesian _polar LIBRARY
  DESTINATION ${PYTHON_PKG_DIR}/geometry
  COMPONENT python)
install(
  FILES cartesian.py  __init__.py  polar.py
  DESTINATION ${PYTHON_PKG_DIR}/geometry
  COMPONENT python)
