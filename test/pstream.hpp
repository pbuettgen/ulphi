// -*- coding: utf-8 -*-
/* 
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and library for
 * computing electromagnetic fields.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

/*!
 *  \file
 *
 *  \brief pipe data to an other program
 *
 *  \author Philipp Büttgenbach
 *
 *  \copyright All rights reserved.
 */


#ifndef PSTREAM_HPP_SuwovEso_
#define PSTREAM_HPP_SuwovEso_

#include <stdio.h>
#include <ios>
#include <iosfwd>

#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/stream.hpp>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

namespace boost_io = boost::iostreams;


//! Sink for writing to an other process
/*!
 *  \tparam Ch is the character type
 *
 *  Have a look at the boost::iostreams library's documentation for
 *  more explanations.
 */
template <typename Ch>
struct OPSink
{
   //! character type
   typedef Ch char_type;

   //! Sink category
   struct category 
      : boost_io::sink_tag,
	boost_io::closable_tag,
	boost_io::flushable_tag {};

#ifdef HAVE_POPEN
   //! popen(3) a given program for writing.
   /*!
    *  \param program is the program which should be popened.
    */
   OPSink(const char* program)
   : file_(popen(program, "w")), ref_count_(new std::size_t)
   {
      *ref_count_ = 0;

      if ( !(this->is_open()) )
	 throw std::ios_base::failure("PSink: popen failed!");
   }

   //! Create an OPSink from another OPSink
   /*!
    *  \param other is the OPSink which shall be copied.
    */
   OPSink(const OPSink<char_type>& other)
      : file_(other.file_), ref_count_(other.ref_count_)
   {
      ++(*ref_count_);
   }

   //! Destruct a OPSink
   /*!
    *  If a opened process is attached to the sink, also pclose(3) the
    *  process.
    */
   ~OPSink() {
      if (! (*ref_count_) ) {
	 this->close();
	 delete ref_count_;
      }
      else {
	 --(*ref_count_);
      }
   }

   //! Determine whether the sink is attached to a process
   bool is_open() const { return nullptr != file_; }

   //! Write something to the process
   /*!
    *  \param s points to the character string which shall be written
    *
    *  \param n is the number of characters to be written
    *
    *  \returns the number of successfully written characters
    */
   std::streamsize write(const char_type* s, std::streamsize n) {
      return std::fwrite(s, sizeof(char_type), n, file_);
   }

   //! Flush the attached stream
   bool flush() { return static_cast<bool>(std::fflush(file_)); }

   //! pclose(3) the attached process
   void close() {
      if (this->is_open())
	 pclose(file_);
      file_=nullptr;
   }

#else
   //! popen(3) a given program for writing.
   /*!
    *  \param program is the program which should be popened.
    */
   OPSink(const *char program) : file_(nullptr), ref_count_(nullptr) {}

   //! Determine whether the sink is attached to a process
   bool is_open() {return false;}

   //! Flush the attached stream
   bool flush() {return false;}

   //! pclose(3) the attached process
   void close() {}

   //! Write something to the process
   /*!
    *  \param s points to the character string which shall be written
    *
    *  \param n is the number of characters to be written
    *
    *  \returns the number of successfully written characters
    */
   std::streamsize write(const char_type* s, std::streamsize n) {return 0;}
#endif

private:
   OPSink<char_type> operator=(const OPSink<char_type>&);

private:
   std::FILE* file_;
   std::size_t* ref_count_;
};

//! Narrow character pipe-stream
typedef boost_io::stream<OPSink<char> > OPStream;

//! Wide character pipe-stream
typedef boost_io::stream<OPSink<wchar_t> > WOPStream;

#endif // PSTREAM_H_
