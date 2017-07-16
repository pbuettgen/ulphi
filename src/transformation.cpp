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

#include <Eigen/Core>

#include "transformation.hpp"

namespace PACKAGE_NS {
   Transformation::~Transformation()
   {
   }

   Shift::~Shift()
   {
   }

   void Shift::Apply(Point& point) const
   {
      using std::get;

      get<0>(point) += get<0>(shift_vector_);
      get<1>(point) += get<1>(shift_vector_);
   }

   std::unique_ptr<Transformation> Shift::Clone() const
   {
      return std::unique_ptr<Transformation>(new Shift(*this));
   }

   Rotate::~Rotate()
   {
   }

   void Rotate::Apply(Point& point) const
   {
      using std::get;

      auto new_point = std::make_tuple(
         get<0>(point) * cos_alpha_ - get<1>(point) * sin_alpha_,
         get<0>(point) * sin_alpha_ + get<1>(point) * cos_alpha_);

      std::swap(point, new_point);
   }

   std::unique_ptr<Transformation> Rotate::Clone() const
   {
      return std::unique_ptr<Transformation>(new Rotate(*this));
   }
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
