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

include(UseMarkdown)

add_markdown_target(
  "white-paper"
  ${CMAKE_CURRENT_BINARY_DIR}
  white-paper.md)

install_markdown_target(
  "white-paper"
  DESTINATION ${CMAKE_INSTALL_DOCDIR}
  COMPONENT "doc")
