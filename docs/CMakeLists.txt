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

set(LATEX_COMPILER_FLAGS "-interaction=batchmode -file-line-error -shell-escape"
  CACHE STRING "Flags passed to latex.")
include(UseLATEX)

# test whether ulphi is installed
include(FindPythonModule)
find_python_module(${PROJECT_NAME})

macro(makeLink src dest)
  add_custom_command(
    OUTPUT ${dest}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E create_symlink ${src} ${dest}
    DEPENDS  ${src}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "linking ${src} to ${dest}")
endmacro()

# Links for pdf includes
macro(MakePdfLinks target)
  set(AllPdfs "")
  foreach(pdf_tex_ ${input_pdf_tex_files})
    string(REGEX REPLACE "\.pdf_tex$" ".pdf" pdf_ ${pdf_tex_})
    get_filename_component(pdf_basename ${pdf_} NAME)
    makeLink(
      ${CMAKE_CURRENT_BINARY_DIR}/${pdf_}
      ${CMAKE_CURRENT_BINARY_DIR}/${pdf_basename})
    list(APPEND AllPdfs ${pdf_basename})
  endforeach()
  add_custom_target(LinkPdfs_${target} DEPENDS ${AllPdfs})
  add_dependencies(LinkPdfs_${target} SvgAsPdf)
endmacro()

add_subdirectory(datasets)
add_subdirectory(images)
add_subdirectory(include)
add_subdirectory(manual)
add_subdirectory(poster)
add_subdirectory(presentation)
add_subdirectory(white-paper)
