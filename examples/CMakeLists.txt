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

set(exampledir ${CMAKE_INSTALL_DOCDIR}/example)

set(example_PYTHON
  coaxial_line.py
  E_I_transformer_core.py
  simple_iron_core.py
  toroidal_wound_core.py
  wound_iron_core.py)

install(
  FILES ${example_PYTHON}
  DESTINATION ${exampledir}
  COMPONENT doc)
