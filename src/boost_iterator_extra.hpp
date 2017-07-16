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
 *  \brief Some extra stuff for boost::iterator
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef BOOST_ITERATOR_EXTRA_HPP_Gloyrowy_
#define BOOST_ITERATOR_EXTRA_HPP_Gloyrowy_

#include <type_traits>

#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace std {
   template <class UnaryFunc, class Iterator, class Reference, class Value>
   struct iterator_traits<boost::transform_iterator<
                             UnaryFunc, Iterator, Reference, Value> >
   {
   private:
      using TransformIterator =
         boost::transform_iterator<UnaryFunc, Iterator, Reference, Value>;

   public:
      using reference =
         typename boost::iterators::detail::ia_dflt_help<
      Reference
      , result_of<const UnaryFunc(typename std::iterator_traits<Iterator>::reference)>
      >::type;

      using value_type =
         typename boost::iterators::detail::ia_dflt_help<
         Value
         , remove_reference<reference>
         >::type;
   };
}

#endif // BOOST_ITERATOR_EXTRA_HPP_Gloyrowy_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
