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

set(Document poster)

set(MainTeXFile ${Document}.tex)
set(ExtraTeXFiles ../include/preamble-common.tex)

foreach(scan_type
    addplot_csv includegraphics include_tex input_pdf_tex input_tex)
  execute_process(
    COMMAND ${PROJECT_SOURCE_DIR}/bin/scan_texfile_includes.py -q --${scan_type} ${MainTeXFile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE ${scan_type}_files)
endforeach()
list(REMOVE_ITEM includegraphics_files "../images/Ulphi-logo.pdf")

MakePdfLinks(${Document})

add_latex_document(${MainTeXFile}
  INPUTS ${ExtraTeXFiles} 
  BIBFILES ../include/literature.bib
  IMAGES ${includegraphics_files}
  DEPENDS DataSets LinkPdfs_${Document}
  USE_BIBLATEX
  )

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${Document}.pdf
  DESTINATION ${CMAKE_INSTALL_DOCDIR}
  COMPONENT doc)
