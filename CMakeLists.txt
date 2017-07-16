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

cmake_minimum_required (VERSION 3.7)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-pipe")
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

project (ulphi VERSION 1.1.0 LANGUAGES CXX C)

set(PROJECT_LINKER_FLAGS "-Wl,--demangle")
set(CMAKE_SHARED_LINKER_FLAGS ${PROJECT_LINKER_FLAGS})
set(CMAKE_MODULE_LINKER_FLAGS ${PROJECT_LINKER_FLAGS})

option(MAINTAINER_MODE "Enable maintainer mode." OFF)
mark_as_advanced(MAINTAINER_MODE)

include(CheckBoostBug9692)
include(GNUInstallDirs)
include(OptimizeForArchitecture)
include(EnableExtraCompilerWarnings)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if("${CMAKE_BUILD_TYPE}" STREQUAL "None")
  set(CMAKE_BUILD_TYPE "Release")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  OptimizeForArchitecture()
endif()

# Checks for libraries.

find_package(PythonInterp 3.2 REQUIRED)
find_package(PythonLibs  3.2 REQUIRED)

## Boost
set(Boost_USE_MULTITHREADED OFF)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(Boost_USE_DEBUG_RUNTIME OFF)
  set(Boost_USE_DEBUG_PYTHON OFF)
else()
  set(Boost_USE_DEBUG_RUNTIME ON)
  set(Boost_USE_DEBUG_PYTHON ON)
endif()
find_package(Boost 1.55 REQUIRED graph log python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
CheckBoostBug9692()

## Eigen
find_package(Eigen3 REQUIRED)

## Loki
find_package(Loki REQUIRED)

# Process subdirectories
add_subdirectory(src)
add_subdirectory(python)
add_subdirectory(test)
if(${MAINTAINER_MODE})
  add_subdirectory(profile)
endif()

# Documentation
add_subdirectory(example)
add_subdirectory(documents)

include(DoxygenTargets)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

add_doxygen(${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  INSTALL_DESTINATION ${CMAKE_INSTALL_DOCDIR}
  INSTALL_COMPONENT doc
  PROJECT_NUMBER ${PROJECT_VERSION}
  NO_PDF)

find_package(ImageMagick COMPONENTS convert REQUIRED)

add_custom_command(
  TARGET doc PRE_BUILD
  COMMAND ${IMAGEMAGICK_CONVERT_EXECUTABLE} -resize 200x55 ${PROJECT_SOURCE_DIR}/documents/images/Ulphi-logo.svg Ulphi-logo.png
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Converting logo from svg to png..."
  VERBATIM)

install(
  FILES debian/${PROJECT_NAME}-doc.doc-base
  DESTINATION ${CMAKE_INSTALL_DATADIR}/doc-base
  RENAME ${PROJECT_NAME}
  COMPONENT doc)

install(
  FILES AUTHORS NEWS README.md
  DESTINATION ${CMAKE_INSTALL_DOCDIR}
  COMPONENT doc)

if(${MAINTAINER_MODE})
  find_program(copyright_update "copyright-update" REQUIRED)
  add_custom_target(
    update-copyright
    COMMAND find \( -name \*.cpp -or -name \*.hpp -or -name Makefile.am
    -or -name COPYING
    -or -name \*.py -or -name \*.tex -or -name \*.bib \)
    -and -not -wholename \*doxygen-doc\* -print0
    | xargs -0 ${copyright_update}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

set(CPACK_COMPONENTS_ALL dev lib python doc)
set(CPACK_GENERATOR DEB)
set(CPACK_INSTALL_CMAKE_PROJECTS ${CMAKE_CURRENT_BINARY_DIR} ${PROJECT_NAME} ALL /)
set(CPACK_MONOLITHIC_INSTALL OFF)
set(CPACK_PACKAGE_CONTACT "Philipp Büttgenbach <Philipp.Buettgenbach@gmx.de>")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "${PROJECT_NAME} -- a library for computing electromagnetic fields")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING.MPL2")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
set(CPACK_RPM_COMPONENT_INSTALL ON)
set(CPACK_STRIP_FILES ON)

# Debian Packages
get_filename_component(
  PYTHON_EXECUTABLE_BASENAME ${PYTHON_EXECUTABLE} NAME)

set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://github.com/pbuettgen/${PROJECT_NAME}")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_DEB_COMPONENT_INSTALL ON)

set(CPACK_DEBIAN_DEV_PACKAGE_SECTION "contrib/devel")
set(CPACK_DEBIAN_DOC_PACKAGE_ARCHITECTURE "all")
set(CPACK_DEBIAN_DOC_PACKAGE_SECTION "contrib/doc")
set(CPACK_DEBIAN_LIB_PACKAGE_SECTION "contrib/libs")
set(CPACK_DEBIAN_LIB_PACKAGE_NAME "lib${PROJECT_NAME}${PROJECT_VERSION_MAJOR}")
set(CPACK_DEBIAN_PYTHON_PACKAGE_SECTION "contrib/python")
set(CPACK_DEBIAN_PYTHON_PACKAGE_NAME
  "${PYTHON_EXECUTABLE_BASENAME}-${PROJECT_NAME}")
set(CPACK_DEBIAN_PYTHON_PACKAGE_DEPENDS
  "${PYTHON_EXECUTABLE_BASENAME} (>= ${PYTHON_VERSION_STRING}), ${PYTHON_EXECUTABLE_BASENAME}-pint (>= 0.7)")

set(CPACK_SOURCE_GENERATOR TGZ ZIP)

set(CPACK_SOURCE_IGNORE_FILES
  "~$"
  "\\\\.svn"
  "\\\\.gitignore$"
  "\\\\.gitmodules$"
  "^${PROJECT_SOURCE_DIR}/.git/"
  "^${PROJECT_SOURCE_DIR}/build/"
  )

include(CPack)

cpack_add_component(dev
  DISPLAY_NAME "development"
  DESCRIPTION "C++ header files and development libraries")
cpack_add_component(lib
  DISPLAY_NAME "library"
  DESCRIPTION "${PROJECT_NAME} library")
cpack_add_component(python
  DISPLAY_NAME "python"
  DESCRIPTION "python modules")
cpack_add_component(doc
  DISPLAY_NAME "documentation"
  DESCRIPTION "${PROJECT_NAME} documentation")
