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
 *  \brief Assemble and solve a linear system from a grid graph
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-07-16 13:55:52 +0200 (So, 16. Jul 2017) $
 *  \copyright MPL 2.0
 */

#ifndef LINEAR_SOLVE_HPP_3014OHqm_
#define LINEAR_SOLVE_HPP_3014OHqm_

#include <Eigen/SparseLU>

#include "config.h"

#ifdef HAVE_BOOST_LOG
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#endif

#include "exception.hpp"
#include "grid_graph.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Set up the equation system's matrix and right hand side
      /*!
       * @tparam GridGraphPointer is the grid graph pointer type
       *
       * @param grid_graph is the graph describing the problem
       * @param system_matrix is the equation system's matrix
       * @param rhs is the equation system's right hand side
       */
      template<typename GridGraphPointer,
               typename Analysis
               = typename GridGraphPointer::element_type::AnalysisT>
      inline void InitLinearEquationSystem(
         GridGraphPointer grid_graph,
         typename Analysis::SystemMatrix& system_matrix,
         typename Analysis::RhsVector& rhs)
      {
         using AssemblingMatrix = typename Analysis::AssemblingMatrix;

         AssemblingMatrix system_triplets;

         rhs.setZero();

         for (auto v = vertices(*grid_graph); v.first != v.second; ++(v.first))
         {
            (*grid_graph)[*(v.first)]->AddToEquationSystem(*(v.first),
                                                           grid_graph, system_triplets, rhs);
         }

         system_matrix.setFromTriplets(system_triplets.begin(),
                                       system_triplets.end());
      }
   }

   //! Do a linear solving process onto a grid
   /*!
    * Useful if all reluctivities are independent from the magnetic
    * flux density.
    *
    * @tparam GridGraphPointer is the grid graph pointer type
    *
    *  \param grid is some grid object
    */
   template<typename GridGraphPointer>
   void LinearSolve(GridGraphPointer grid_graph)
   {
      using Analysis = typename GridGraphPointer::element_type::AnalysisT;
      using SystemMatrix = typename Analysis::SystemMatrix;
      using RhsVector = typename Analysis::RhsVector;

      const EquationNumber number_of_equations =
         NumberOfEquations(grid_graph);

      SystemMatrix system_matrix(number_of_equations, number_of_equations);
      RhsVector rhs(number_of_equations), solution(number_of_equations);

      detail::InitLinearEquationSystem(grid_graph, system_matrix, rhs);

      Eigen::SparseLU<
         SystemMatrix,
         Eigen::COLAMDOrdering<typename SystemMatrix::Index> > solver(
            system_matrix);
      auto error_message = solver.lastErrorMessage();
      if ("" != error_message)
         throw RuntimeError(error_message);
      solution = solver.solve(rhs);

      SetPotential(grid_graph, solution);
      UpdateMagneticFluxDensity(grid_graph);
   }

   //! Solve a nonlinear system by simple iterations
   /*!
    * The system is resolved until convergence seems to be reached.
    * May oszillate.  In most cases Newton is probably a better choice
    * for solving nonlinear systems.
    *
    * @tparam GridGraphPointer is the grid graph pointer type
    *
    * @param grid_graph pointer to a grid graph.
    *
    * @param epsilon upper limit for a remaining error.
    *
    * @param[in, out] max_iter maximum number of iterations, contains
    * the actual number of iterations on return
    */
   template <typename GridGraphPointer>
   void IterativeSolve(
      GridGraphPointer grid_graph, double epsilon=1e-4,
      std::size_t& max_iter = std::numeric_limits<std::size_t>::max()/2)
   {
      using Analysis = typename GridGraphPointer::element_type::AnalysisT;
      using SystemMatrix = typename Analysis::SystemMatrix;
      using RhsVector = typename Analysis::RhsVector;

      const EquationNumber number_of_equations =
         NumberOfEquations(grid_graph);

      SystemMatrix system_matrix(number_of_equations, number_of_equations);
      RhsVector
         rhs(number_of_equations),
         solution(number_of_equations),
         previous_solution(number_of_equations);

      previous_solution.setZero();

      std::size_t iter_count = max_iter;
      double error = 0.;

#ifdef HAVE_BOOST_LOG
      boost::format iter_message("Iteration: %1% -- error: %2%");
#endif

      do {
         detail::InitLinearEquationSystem(grid_graph, system_matrix, rhs);

         Eigen::SparseLU<SystemMatrix, Eigen::COLAMDOrdering<
            typename SystemMatrix::Index> > solver(system_matrix);
         auto error_message = solver.lastErrorMessage();
         if ("" != error_message)
            throw RuntimeError(error_message);
         solution = solver.solve(rhs);

         error = (solution-previous_solution).norm();

         previous_solution = solution;

         SetPotential(grid_graph, solution);
         UpdateMagneticFluxDensity(grid_graph);

#ifdef HAVE_BOOST_LOG
         BOOST_LOG_TRIVIAL(info) << (iter_message % (max_iter-iter_count)
                                     % error).str();
#endif
      } while(--iter_count and error > epsilon);

      max_iter -= iter_count;
   }
}

#endif // LINEAR_SOLVE_HPP_3014OHqm_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
