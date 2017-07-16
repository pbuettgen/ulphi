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

include(CheckCXXSourceCompiles)

macro(CheckBoostBug9692)
  include_directories(${Boost_INCLUDE_DIR})
  check_cxx_source_compiles("
// Check whether bi-dir-csr-graph supports global graph properties
#include <boost/graph/compressed_sparse_row_graph.hpp>
struct VertexProperty {};
struct EdgeProperty {};
struct GraphProperty {int x;};

int main() {
using namespace boost;
compressed_sparse_row_graph<bidirectionalS, VertexProperty, EdgeProperty, GraphProperty> csr_graph;
csr_graph[graph_bundle].x=1;
}" HAVE_CSR_BIDIR_GRAPH_PROPERTY_ACCESS)
endmacro()
