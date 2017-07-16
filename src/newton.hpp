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
 *  \brief Solve a nonlinear problem with Newton's method.
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-07-16 13:55:52 +0200 (So, 16. Jul 2017) $
 *  \copyright MPL 2.0
 */

#ifndef NEWTON_HPP_Q8Idx7wQ_
#define NEWTON_HPP_Q8Idx7wQ_

#include <cstdlib>
#include <fstream>
#include <limits>

#include <Eigen/SparseLU>

#include "config.h"
#ifdef HAVE_BOOST_LOG
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#endif

#include "analysis.hpp"
#include "exception.hpp"
#include "grid_graph.hpp"
#include "vertex.hpp"
#include "types.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Update the equation system after a newton step
      /*!
       *  \tparam GridGraphPointer is the grid graph pointer type
       *
       *  \param grid_graph is a GridGraph
       *  \param system_matrix is the equation system's matrix
       *  \param jacobian is the jacobi matrix
       *  \param rhs is the equation system's right hand side
       */
      template<typename GridGraphPointer,
               typename Analysis =
               typename GridGraphPointer::element_type::AnalysisT>
      void MatrixUpdate(
         GridGraphPointer grid_graph,
         typename Analysis::SystemMatrix& system_matrix,
         typename Analysis::SystemMatrix& jacobian,
         typename Analysis::RhsVector& rhs)
      {
         using AssemblingMatrix = typename Analysis::AssemblingMatrix;

         rhs.setZero();

         AssemblingMatrix system_triplets, jacobi_triplets;

         for (auto v = vertices(*grid_graph); v.first != v.second;
              ++(v.first)) {
            (*grid_graph)[*(v.first)]->AddToEquationSystem(*(v.first),
                                                           grid_graph, system_triplets, jacobi_triplets, rhs);
         }

         system_matrix.setFromTriplets(system_triplets.begin(),
                                       system_triplets.end());
         jacobian.setFromTriplets(jacobi_triplets.begin(),
                                  jacobi_triplets.end());
      }
   }

   //! Solve a nonlinear equation system using newton's method
   /*!
    *  \tparam GridGraphPointer is the grid graph pointer type
    *
    *  \param grid_graph points to a GridGraph
    *
    *  \param epsilon defines the precision which shall be reached
    *
    *  \param[in,out] max_iter is the maximum number of iterations.
    *  Before the function returns it is set to the actual number of
    *  iterations.
    *
    *  \throw RuntimeError if solving the equation system fails.
    *
    *  \throw RuntimeError if the damping becomes too strong.
    */
   template<typename GridGraphPointer>
   void Newton(
      GridGraphPointer grid_graph, double epsilon = 1e-4,
      std::size_t max_iter = 36)
   {
      using Analysis = typename GridGraphPointer::element_type::AnalysisT;
      using SystemMatrix = typename Analysis::SystemMatrix;
      using RhsVector = typename Analysis::RhsVector;

      CheckValidGridGraphPointer(grid_graph);

      const EquationNumber number_of_equations
         = NumberOfEquations(grid_graph);

      RhsVector rhs(number_of_equations), delta(number_of_equations), result(
         number_of_equations);
      SystemMatrix system_matrix(number_of_equations, number_of_equations),
         jacobian(number_of_equations, number_of_equations);

      result.setZero();

      std::size_t count = max_iter;
      double error = .0;
      double last_error = std::numeric_limits<double>::max();
      double damping_factor = 1.;
      constexpr double min_damping_factor = 1e-4;

#ifdef HAVE_BOOST_LOG
      boost::format iter_message(
         "Newton: iteration: %1% -- norm: %2% -- rel. error: %3% %% -- damping: %4%");
#endif

      do {
         MatrixUpdate(grid_graph, system_matrix, jacobian, rhs);

         rhs -= system_matrix * result;

         Eigen::SparseLU<
            SystemMatrix,
            Eigen::COLAMDOrdering<typename SystemMatrix::Index> > solver(
               jacobian);
         auto error_message = solver.lastErrorMessage();
         if ("" != error_message)
            throw RuntimeError(error_message);
         delta = solver.solve(rhs);

         error = delta.norm();

         while (damping_factor * error >= last_error) {
            damping_factor /= 2.;
            if (min_damping_factor > damping_factor) {
               throw RuntimeError(
                  "Newton: Damping too strong. No convergence!");
            }
         }

         result += damping_factor * delta;

         last_error = error;

         SetPotential(grid_graph, result);
         UpdateMagneticFluxDensity(grid_graph);

#ifdef HAVE_BOOST_LOG
         double rel_error = 100.*error/result.cwiseAbs().maxCoeff();
         BOOST_LOG_TRIVIAL(info) << (iter_message % (max_iter-count) % error
                                     % rel_error % damping_factor).str();
#endif
         damping_factor = std::min(1., 2. * damping_factor);
      } while (--count and error > epsilon);

      // max_iter -= count;
   }
}

#endif // NEWTON_HPP_Q8Idx7wQ_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
