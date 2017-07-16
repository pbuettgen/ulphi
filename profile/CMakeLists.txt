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

find_program(gprof gprof)
if(x${gprof} STREQUAL x  )
  message(ERROR "Please install gprof for profiling")
endif()

# enable profiling
foreach(profile_switch -pg)
  check_cxx_compiler_flag(${profile_switch} cxx_has_this_switch)
  if(${cxx_has_this_switch})
    add_compile_options(${profile_switch})
    break()
  endif()
endforeach()

add_executable(
  E_I_transformer_core EXCLUDE_FROM_ALL
  E-I-transformer-core.cpp
  ${PROJECT_SOURCE_DIR}/src/constant_reluctivity.cpp
  ${PROJECT_SOURCE_DIR}/src/linear_vertex_distribution.cpp
  ${PROJECT_SOURCE_DIR}/src/scalar_reluctivity.cpp
  ${PROJECT_SOURCE_DIR}/src/spline_reluctivity.cpp
  ${PROJECT_SOURCE_DIR}/src/vertex_distribution.cpp)

target_compile_options(
  E_I_transformer_core PRIVATE -include magnetostatic_analysis.hpp)
target_compile_definitions(
  E_I_transformer_core PRIVATE INSTANTIATE_)

target_link_libraries(
  E_I_transformer_core ${LOKI_LIBRARY} ${BOOST_LOG_LIBRARY})

add_custom_command(
  OUTPUT gmon.out
  COMMAND ./E_I_transformer_core
  DEPENDS E_I_transformer_core
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(
  OUTPUT report.txt
  COMMAND ${gprof} E_I_transformer_core gmon.out > report.txt
  DEPENDS gmon.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(profile DEPENDS report.txt)
