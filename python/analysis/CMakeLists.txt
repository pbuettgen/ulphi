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

foreach(analysis_type static harmonic)
  PYTHON_ADD_MODULE(
    _magneto${analysis_type} analysis.cpp
    $<TARGET_OBJECTS:lib${PROJECT_NAME}_base>
    $<TARGET_OBJECTS:lib${PROJECT_NAME}_analysis_magneto${analysis_type}>
    $<TARGET_OBJECTS:lib${PROJECT_NAME}_materials_basic>
    $<TARGET_OBJECTS:lib${PROJECT_NAME}_transformations>)
  target_compile_definitions(
    _magneto${analysis_type} PRIVATE MODULE_NAME=_magneto${analysis_type})
  target_compile_options(
    _magneto${analysis_type} PRIVATE
    -include magneto${analysis_type}_analysis.hpp)
  target_link_libraries(_magneto${analysis_type} ${LOKI_LIBRARY})
  target_link_libraries(
    _magneto${analysis_type} ${Boost_LOG_LIBRARY})
  install(
    TARGETS _magneto${analysis_type} LIBRARY
    DESTINATION ${PYTHON_PKG_DIR}/analysis
    COMPONENT python)
  install(
    FILES magneto${analysis_type}.py
    DESTINATION ${PYTHON_PKG_DIR}/analysis
    COMPONENT python)
endforeach()

install(
  FILES __init__.py
  DESTINATION ${PYTHON_PKG_DIR}/analysis
  COMPONENT python)
