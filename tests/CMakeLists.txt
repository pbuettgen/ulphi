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

find_path(
  GTEST_SRC_DIR
  NAMES CMakeLists.txt
  HINTS /usr/src/googletest/googletest /usr/src/gtest
  DOC "Google test library's source directory"
  NO_DEFAULT_PATH
  REQUIRED)

add_subdirectory(${GTEST_SRC_DIR}  gtest)

find_package(GTest)

enable_testing()

set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

include_directories(
  ${PROJECT_BINARY_DIR}/src
  ${PROJECT_SOURCE_DIR}/src
  ${EIGEN3_INCLUDE_DIR}
  ${GTEST_INCLUDE_DIRS})

link_directories(
  ${CMAKE_CURRENT_BINARY_DIR}/gtest)

add_executable(
  checker
  ${CMAKE_CURRENT_SOURCE_DIR}/checker.cpp)
add_dependencies(
  checker
  ${PROJECT_NAME} gtest_main Threads::Threads)
target_link_libraries(
  checker
  gtest_main ${PROJECT_NAME} ${LOKI_LIBRARY})

GTEST_ADD_TESTS(checker "" "AUTO")

add_executable(
  checker_python
  ${CMAKE_CURRENT_SOURCE_DIR}/checker_python.cpp)
add_dependencies(
  checker_python
  ${PROJECT_NAME} gtest_main Threads::Threads)
target_include_directories(
  checker_python PRIVATE ${PROJECT_SOURCE_DIR} ${PYTHON_INCLUDE_DIR})
target_link_libraries(
  checker_python
  gtest_main
  ${Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_LIBRARY_RELEASE}
  ${PYTHON_LIBRARY})

GTEST_ADD_TESTS(checker_python "" "AUTO")
