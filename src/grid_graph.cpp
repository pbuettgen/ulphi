// -*- coding: utf-8 -*-
/*
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and library for computing
 * electromagnetic fields.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "grid_graph-fwd.hpp"
#include "grid_graph.hpp"
#include "standard_vertex.hpp"
#include "neumann_vertex.hpp"
#include "dirichlet_vertex.hpp"

namespace PACKAGE_NS {
   namespace fi {
      //template struct GridGraph<TimeTransientAnalysis<> >;
      //template struct GridGraph<TimeHarmonicAnalysis<> >;
   }

#define VERTEX_TEMPLATES                                \
   template struct StandardVertex<AnalysisType>;        \
   template struct NeumannVertex<AnalysisType>;         \
   template struct ConstDirichletVertex<AnalysisType>;	\
   template struct VarDirichletVertex<AnalysisType>;

   namespace fi {
      //VERTEX_TEMPLATES
   }
}

#ifdef VERTEX_TEMPLATES
#undef VERTEX_TEMPLATES
#endif

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
