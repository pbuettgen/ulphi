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

#
# Inkscape graphics
#
find_program(INKSCAPE_EXECUTABLE inkscape)

file(GLOB svgfiles_ RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.svg)

set(all_pdf_from_svg "")
set(all_pdf_tex "")

foreach(svgf_ ${svgfiles_})
  string(REGEX REPLACE "\.svg$" ".pdf" pdf_ ${svgf_})
  string(REGEX REPLACE "\.svg$" ".pdf_tex" pdf_tex_ ${svgf_})
  list(APPEND all_pdf_from_svg ${pdf_})
  list(APPEND all_pdf_tex ${pdf_tex_})
  add_custom_command(
    OUTPUT ${pdf_} ${pdf_tex_}
    COMMAND ${INKSCAPE_EXECUTABLE} -z -D --export-latex --export-pdf=${pdf_} ${CMAKE_CURRENT_SOURCE_DIR}/${svgf_}
    DEPENDS ${svgf_}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Exporting ${svgf_} ...")
endforeach(svgf_)

add_custom_target(SvgAsPdf DEPENDS ${all_pdf_tex} ${all_pdf_from_svg})

#
# Sketches
#
find_program(M4_EXECUTABLE m4)
find_path(M4_INCLUDE_DIR
  NAMES "for.m4" "foreach.m4"
  PATHS "/usr/share/doc/m4/examples")
find_program(SKETCH_EXECUTABLE sketch)

file(GLOB skfiles_ RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.sk)

set(all_tex_from_sk "")
foreach(skf_ ${skfiles_})
    string(REGEX REPLACE "\.sk$" ".tex" tex_ ${skf_})
    list(APPEND all_tex_from_sk ${tex_})
    add_custom_command(
      OUTPUT ${tex_}
      COMMAND ${M4_EXECUTABLE} -I ${M4_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${skf_} | ${SKETCH_EXECUTABLE} -o ${tex_}
      DEPENDS ${skf_}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Processing ${skf_} ...")
endforeach()

add_custom_target(SkAsTex DEPENDS ${all_tex_from_sk})

#
# Other images
#
file(
  GLOB other_images_
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.png *.jpg *.pdf)

foreach(if_ ${other_images_})
  makeLink(
    ${CMAKE_CURRENT_SOURCE_DIR}/${if_}
    ${CMAKE_CURRENT_BINARY_DIR}/${if_})
endforeach(if_)
