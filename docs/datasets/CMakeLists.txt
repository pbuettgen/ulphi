# -*- coding: utf-8 -*-
#
# Copyright © 2012-2016 Philipp Büttgenbach
#
# This file is part of ulphi, a CAE tool and library for
# computing electromagnetic fields.
#
# All rights reserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#

string(TOUPPER ${PROJECT_NAME} project_name_upper)

file(GLOB datafiles1 RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "M??0-*.csv")
foreach(df1 ${datafiles1})
  makeLink(
    ${CMAKE_CURRENT_SOURCE_DIR}/${df1}
    ${CMAKE_CURRENT_BINARY_DIR}/${df1})
endforeach()

if(${PY_${project_name_upper}})
  foreach(df3 coaxial_line simple_iron_core)
    add_custom_command(
      OUTPUT ${df3}.csv
      COMMAND ${PROJECT_SOURCE_DIR}/example/${df3}.py -f ${CMAKE_CURRENT_BINARY_DIR}/{df3}.csv
      DEPENDS ${PROJECT_SOURCE_DIR}/example/${df3}.py)
  endforeach()
  foreach(df4 E_I_transformer_core wound_iron_core)
    add_custom_command(
      OUTPUT ${df4}_freq.csv ${df4}_time.csv
      COMMAND ${PROJECT_SOURCE_DIR}/example/${df4}.py -f ${CMAKE_CURRENT_BINARY_DIR}/${df4}
      DEPENDS ${PROJECT_SOURCE_DIR}/example/${df4}.py)
  endforeach()
else()
  foreach(df2
      coaxial_line.csv
      E_I_transformer_core_freq.csv
      E_I_transformer_core_time.csv
      simple_iron_core.csv
      wound_iron_core_freq.csv
      wound_iron_core_time.csv)
    makeLink(
      ${CMAKE_CURRENT_SOURCE_DIR}/${df2}
      ${CMAKE_CURRENT_BINARY_DIR}/${df2})
  endforeach()
endif()

if(${PY_${project_name_upper}})
  foreach(type arctan sqrt spline)
    add_custom_command(
      OUTPUT M400_50A-${type}.csv
      COMMAND ${PROJECT_SOURCE_DIR}/bin/reluctivity-models.py -r ${type}
      DEPENDS ${PROJECT_SOURCE_DIR}/bin/reluctivity-models.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endforeach()
else()
  foreach(type arctan sqrt spline)
    makeLink(
      ${CMAKE_CURRENT_SOURCE_DIR}/M400_50A-${type}.csv
      ${CMAKE_CURRENT_BINARY_DIR}/M400_50A-${type}.csv)
  endforeach()
endif()

add_custom_target(
  DataSets
  DEPENDS ${datafiles1} coaxial_line.csv
  E_I_transformer_core_freq.csv E_I_transformer_core_time.csv
  simple_iron_core.csv
  wound_iron_core_freq.csv wound_iron_core_time.csv
  M400_50A-arctan.csv
  M400_50A-sqrt.csv
  M400_50A-spline.csv)
