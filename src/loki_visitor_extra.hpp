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

/*!
 *  \file
 *
 *  \brief extras for the loki visitor library
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef LOKI_VISITOR_EXTRA_HPP_AGXljNV4_
#define LOKI_VISITOR_EXTRA_HPP_AGXljNV4_

#include "config.h"

#ifndef HAVE_LIBLOKI
#error "The loki library is now required"
#endif

#include <loki/Visitor.h>

/*!
 * \def LOKI_DEFINE_VISITABLE()
 * \ingroup VisitorGroup
 *
 * Put it in every class that you want to make visitable
 * (in addition to deriving it from BaseVisitable<R>)
 */
#ifdef LOKI_DEFINE_VISITABLE
#undef LOKI_DEFINE_VISITABLE
#define LOKI_DEFINE_VISITABLE()					\
    virtual ReturnType Accept(::Loki::BaseVisitor& guest)	\
    { return this->AcceptImpl(*this, guest); }
#endif

/*!
 * \def LOKI_DEFINE_CONST_VISITABLE()
 * \ingroup VisitorGroup
 *
 * Put it in every class that you want to make visitable by const
 * member functions (in addition to deriving it from BaseVisitable<R>)
 */
#ifdef LOKI_DEFINE_CONST_VISITABLE
#undef LOKI_DEFINE_CONST_VISITABLE
#define LOKI_DEFINE_CONST_VISITABLE()				\
    virtual ReturnType Accept(::Loki::BaseVisitor& guest) const \
    { return this->AcceptImpl(*this, guest); }
#endif

#endif // LOKI_VISITOR_EXTRA_HPP_AGXljNV4_
