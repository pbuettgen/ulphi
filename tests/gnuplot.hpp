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
 *  \brief Stream data to gnuplot for plotting
 *
 *  \author Philipp Büttgenbach
 *
 *  \copyright All rights reserved.
 */


#ifndef GNUPLOT_HPP_FurvAkci_
#define GNUPLOT_HPP_FurvAkci_

#include <string>
#include <fstream>

#include <boost/scope_exit.hpp>

#include "magnetic_material.hpp"
#include "pstream.hpp"


//! Display a material's reluctivity using gnuplot
/*!
 *  \tparam BIterator is an iterator for a sequence of flux density
 *  values
 *
 *  \param bbegin points to the flux density sequence's start.
 *
 *  \param bend points to the flux density sequence's end.
 *
 *  \param material is the material whose reluctivity should be ploted.
 *
 *  \note This function is only available when popen is available on
 *  your system.
 */
template <typename BIterator>
void GnuplotReluctivity(
   const BIterator& bbegin, const BIterator& bend,
   const IsotropicMaterial& material)
{
   std::string datafilename(std::tmpnam(nullptr));

   {
      std::ofstream datafile(datafilename.c_str());

      for (BIterator b=bbegin; bend!=b; ++b)
      {
	 datafile 
	    << boost::units::quantity_cast<double>(*b) << "\t"
	    << boost::units::quantity_cast<double>(material.Reluctivity(*b))
	    << "\n";
      }

      datafile.close();
   }

   BOOST_SCOPE_EXIT( (&datafilename) )
   {
      std::remove(datafilename.c_str());
   } BOOST_SCOPE_EXIT_END


   OPStream gnuplot("gnuplot");

   gnuplot << "set grid\nplot \"" << datafilename 
	   << "\" with linespoints\npause mouse" << std::endl;
}


//! Display a material's derived reluctivity
/*!
 *  Display a material's magnetic reluctivity derived with respect to
 *  the magnetic field density B using gnuplot.
 *
 *  \tparam BIterator is an iterator for a sequence of flux density
 *  values
 *
 *  \param bbegin points to the flux density sequence's start.
 *
 *  \param bend points to the flux density sequence's end.
 *
 *  \param material is the material whose reluctivity should be ploted.
 *
 *  \note This function is only available when popen is available on
 *  your system.
 */
template <typename BIterator>
void GnuplotReluctivityDerived(
   const BIterator& bbegin, const BIterator& bend,
   const IsotropicMaterial& material)
{
   std::string datafilename(std::tmpnam(nullptr));

   {
      std::ofstream datafile(datafilename.c_str());

      for (BIterator b=bbegin; bend!=b; ++b)
      {
	 datafile 
	    << boost::units::quantity_cast<double>(*b) << "\t"
	    << boost::units::quantity_cast<double>(material.ReluctivityDerived(*b))
	    << "\n";
      }

      datafile.close();
   }
   FileRemover fr(datafilename);

   OPStream gnuplot("gnuplot");

   gnuplot << "set grid\nplot \"" << datafilename 
	   << "\" with linespoints\npause mouse" << std::endl;
}


#endif // GNUPLOT_H_
