// -*- coding: utf-8 -*-
/* 
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and
 * library for computing electromagnetic fields.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOOST_PYTHON_EXTRA_HPP_GruGridA_
#define BOOST_PYTHON_EXTRA_HPP_GruGridA_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/string.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/preprocessor/stringize.hpp>


namespace PACKAGE_NS {
   namespace python {
      namespace bp = boost::python;

      template <typename Type>
      using BPClassSharedPtr = bp::class_<
         Type, std::shared_ptr<Type>, boost::noncopyable>;
   }
}

#endif
