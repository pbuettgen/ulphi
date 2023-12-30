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

macro(python_add_boost_module module_name)
  PYTHON_ADD_MODULE(_${module_name} ${module_name}.cpp)
  target_compile_definitions(
    _${module_name} PRIVATE MODULE_NAME=_${module_name})
endmacro()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}/src
  ${PROJECT_SOURCE_DIR}/src
  ${Boost_INCLUDE_DIR}
  ${EIGEN3_INCLUDE_DIR}
  ${LOKI_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIR})

link_libraries(
  ${Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_LIBRARY_RELEASE}
  ${PYTHON_LIBRARY})

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(prefix='', plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

file(TO_CMAKE_PATH "${PYTHON_SITE_DIR}" PYTHON_SITE_DIR)

set(PYTHON_PKG_DIR ${PYTHON_SITE_DIR}/${PROJECT_NAME}
  CACHE STRING "python package dir")

add_subdirectory(analysis)
add_subdirectory(geometry)

python_add_boost_module(type_conversions)

PYTHON_ADD_MODULE(
  _materials materials.cpp
  $<TARGET_OBJECTS:lib${PROJECT_NAME}_materials_basic>
  $<TARGET_OBJECTS:lib${PROJECT_NAME}_materials_models>
  $<TARGET_OBJECTS:lib${PROJECT_NAME}_materials_extra>)
target_compile_definitions(_materials PRIVATE MODULE_NAME=_materials)
target_compile_options(_materials
  PRIVATE "-include${PROJECT_SOURCE_DIR}/python/boost/parameter/python.hpp")

PYTHON_ADD_MODULE(
  vertex_distributions vertex_distributions.cpp
  $<TARGET_OBJECTS:lib${PROJECT_NAME}_vertex_distributions>)
target_compile_definitions(
  vertex_distributions PRIVATE MODULE_NAME=vertex_distributions)

install(
  TARGETS vertex_distributions LIBRARY
  DESTINATION ${PYTHON_PKG_DIR}
  COMPONENT python)

install(
  TARGETS _type_conversions LIBRARY
  DESTINATION ${PYTHON_PKG_DIR}
  COMPONENT python)
install(
  FILES __init__.py _quantity.py
  DESTINATION ${PYTHON_PKG_DIR}
  COMPONENT python)

install(
  TARGETS _materials LIBRARY
  DESTINATION ${PYTHON_PKG_DIR}
  COMPONENT python)
install(
  FILES materials.py
  DESTINATION ${PYTHON_PKG_DIR}
  COMPONENT python)
