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
 *  \brief Exception objects
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef EXCEPTION_HPP_6F6JTwQm_
#define EXCEPTION_HPP_6F6JTwQm_

#include <stdexcept>

#include "config.h"

namespace PACKAGE_NS {

   //! Exception class to report arguments outside of expected range
   /*!
    *  neumann uses its own class in order to distinguish neumann
    *  errors from standard library errors.  In all other aspects this
    *  class behaves exactly like the standard one.
    */
   struct OutOfRange: std::out_of_range {
      //! Initialize a OutOfRange exception
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit OutOfRange(const std::string& what_arg)
         : std::out_of_range(what_arg)
      {
      }

      //! Initialize a OutOfRange exception
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit OutOfRange(const char* what_arg)
         : std::out_of_range(what_arg)
      {
      }
   };

   //! Exception class to indicate conditions only detectable at run time
   /*!
    *  neumann uses its own class in order to distinguish neumann
    *  errors from standard library errors.  In all other aspects this
    *  class behaves exactly like the standard one.
    */
   struct RuntimeError: std::runtime_error {
      //! Initialization
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit RuntimeError(const std::string& what_arg)
         : std::runtime_error(what_arg)
      {
      }

      //! Initialization
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit RuntimeError(const char* what_arg)
         : std::runtime_error(what_arg)
      {
      }
   };

   //! exception class to report invalid arguments
   struct InvalidArgument: std::invalid_argument {
      //! Initialization
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit InvalidArgument(const std::string& what_arg)
         : std::invalid_argument(what_arg)
      {
      }

      //! Initialization
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit InvalidArgument(const char* what_arg)
         : std::invalid_argument(what_arg)
      {
      }
   };

   //! exception class to report domain errors
   struct DomainError: std::domain_error {
      //! Initialization
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit DomainError(const std::string& what_arg)
         : std::domain_error(what_arg)
      {
      }

      //! Initialization
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit DomainError(const char* what_arg)
         : std::domain_error(what_arg)
      {
      }
   };

   //! Logic error
   struct LogicError : std::logic_error {
      //! Initialization
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit LogicError(const std::string& what_arg)
         : std::logic_error(what_arg)
      {
      }

      //! Initialization
      /*!
       * @param what_arg gives a description of the error occurred
       */
      explicit LogicError(const char* what_arg)
         : std::logic_error(what_arg)
      {
      }
   };
}

#endif // EXCEPTION_HPP_6F6JTwQm_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
