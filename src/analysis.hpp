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
 *  \brief analysis types
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef ANALYSIS_HPP_g8jiDA3Z_
#define ANALYSIS_HPP_g8jiDA3Z_

#include <complex>
#include <vector>

#include "config.h"

#include <boost/math/constants/constants.hpp>
#include <boost/units/systems/si.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "boost_units_extra.hpp"
#include "exception.hpp"
#include "types.hpp"

namespace PACKAGE_NS {
   namespace detail {
      //! Data types common to time-domain and frequency-domain analysis
      /*!
       * \tparam DType is something like float or double for a
       * time-domain analysis and std::complex for a frequency-domain
       * analysis.
       */
      template<typename DType>
      struct CommonAnalysisTypes {
         //! Data type used for electromagnetic quantities
         /*!
          *  This data type is used for processing all magnetic quantities
          *  expect the material properties which are always real.  For a
          *  time transient analysis this should be some real floating point
          *  data type.
          */
         using DataType = DType;

         //! Data type for magnetic vector potentials
         using MagneticVectorPotential = boost::units::quantity<
            boost::units::si::magnetic_vector_potential, DataType>;

         //! Data type for current densities
         using CurrentDensity = boost::units::quantity<
            boost::units::si::current_density, DataType>;

         //! Data type for current
         using Current = boost::units::quantity<
            boost::units::si::current, DataType>;

         //! Data type for magnetic flux densities
         using MagneticFluxDensity = boost::units::quantity<
            boost::units::si::magnetic_flux_density, DataType>;

         //! Magnetic field intensity
         using MagneticFieldIntensity = boost::units::quantity<
            boost::units::si::magnetic_field_intensity, DataType>;

         //! Magnetic voltage
         using MagneticVoltage = Current;

         //! The equation system's right hand side vector type
         using RhsVector = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;

         //! Matrix type used for matrix assembling
         using AssemblingMatrix = std::vector<Eigen::Triplet<DataType> >;

         //! Matrix type used for system solving
         using SystemMatrix = Eigen::SparseMatrix<DataType, Eigen::ColMajor>;
      };
   }

   //! Datatypes and properties for a time-domain analysis
   /*!
    *  Definition of all properties and data types specific to a
    *  magnetostatic or time transient analysis.
    */
   template<typename Float = double>
   struct TimeTransientAnalysis: detail::CommonAnalysisTypes<Float> {
   private:
      using CommonTypes = detail::CommonAnalysisTypes<Float>;

   public:
      //! Floating point type
      using Real = Float;
      using typename CommonTypes::AssemblingMatrix;
      using typename CommonTypes::CurrentDensity;
      using typename CommonTypes::DataType;
      using typename CommonTypes::MagneticFluxDensity;
      using typename CommonTypes::MagneticVectorPotential;
      using typename CommonTypes::RhsVector;
      using typename CommonTypes::SystemMatrix;

      //! Time data type
      using Time = boost::units::quantity<boost::units::si::time, DataType>;

      //! Number of previously computed potential values to keep
      static constexpr unsigned kTimeIntegratorDegree = 0;
   };

   //! Datatypes and properties for a frequency domain analysis
   /*!
    *  Definition of all properties and data types specific to a time
    *  harmonic analysis.
    */
   template<typename Float = double>
   struct MagnetoHarmonicAnalysis
      : detail::CommonAnalysisTypes<std::complex<Float> > {
   private:
      using CommonTypes = detail::CommonAnalysisTypes<std::complex<Float> >;

   public:
      //! Floating point type
      using Real = Float;
      using typename CommonTypes::AssemblingMatrix;
      using typename CommonTypes::CurrentDensity;
      using typename CommonTypes::DataType;
      using typename CommonTypes::MagneticFluxDensity;
      using typename CommonTypes::MagneticVectorPotential;
      using typename CommonTypes::RhsVector;
      using typename CommonTypes::SystemMatrix;

      //! Number of previously computed potential values to keep
      static constexpr unsigned kTimeIntegratorDegree = 0;

      static constexpr auto kID = "MagnetoHarmonicAnalysis";

      //! Time data type
      /*!
       * In the case of a time harmonic analysis this is actually the
       * frequency
       */
      using Time = boost::units::quantity<boost::units::si::frequency, double>;
   };
}

#endif // ANALYSIS_HPP_g8jiDA3Z

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
