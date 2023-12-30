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

file(GLOB AllFrames RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "frame??.tex")

set(ExtraTeXFiles ${AllFrames} frames.tex ../include/preamble-common.tex)

foreach(scan_type
    addplot_csv includegraphics include_tex input_pdf_tex input_tex)
  execute_process(
    COMMAND ${PROJECT_SOURCE_DIR}/bin/scan_texfile_includes.py -q --${scan_type} frames.tex
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE ${scan_type}_files)
endforeach()
list(REMOVE_ITEM includegraphics_files "../images/Ulphi-logo.pdf")

MakePdfLinks(presentation)

foreach(Document presentation handout)

  set(MainTeXFile ${Document}.tex)

  add_latex_document(
    ${MainTeXFile}
    INPUTS ${ExtraTeXFiles}
    IMAGES ${includegraphics_files}
    DEPENDS LinkPdfs_presentation DataSets)

  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${Document}.pdf
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
    COMPONENT doc)

endforeach()
