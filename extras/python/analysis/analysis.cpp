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

#include <memory>
#include <utility>
#include <iostream>

#include "config.h"

#include <boost/iterator/transform_iterator.hpp>

#include "boost_python_extra.hpp"

#include "dirichlet_condition.hpp"
#include "exp_vertex_distribution.hpp"
#include "grid_assembly.hpp"
#include "grid_compiler.hpp"
#include "grid_patch.hpp"
#include "linear_vertex_distribution.hpp"
#include "linear_solve.hpp"
#include "neumann_condition.hpp"
#include "newton.hpp"


namespace PACKAGE_NS {
   namespace python {
      // -----------------
      // *** GridPatch ***
      // -----------------

      template <typename CoordinateSystem, typename Analysis>
      struct GridPatchWrap : GridPatch<CoordinateSystem, Analysis>
      {
         using Base = GridPatch<CoordinateSystem, Analysis>;
         using Time = typename Analysis::Time;
         using CurrentDensity = typename Analysis::CurrentDensity;
         using Current = typename Analysis::Current;
         using ShapePointer = ConstShapePointer<CoordinateSystem>;

         using Base::Base;
         using Base::SetLoad;

         void SetLoad(ShapePointer shape, bp::object load) {
            bp::extract<CurrentDensity> current_density_extractor(load);
            if (current_density_extractor.check())
               return SetLoad(shape, current_density_extractor());

            bp::extract<Current> current_extractor(load);
            if (current_extractor.check())
               return SetLoad(shape, current_extractor());

            // probably a function of time
            auto test_call_result = load(Time::from_value(0.));

            bp::extract<Current> check_test_call_result(test_call_result);
            if (check_test_call_result.check())
               // A function of time returning some current
               return this->SetLoad(shape, [load](Time t) -> Current {
                     return bp::extract<Current>(load(t));
                  });

            // Last chance: A function of time returning some current
            // density
            return this->SetLoad(shape, [load](Time t) -> CurrentDensity {
                  return bp::extract<CurrentDensity>(load(t));
               });
         }
      };

      template<typename CoordinateSystem, typename Analysis>
      std::shared_ptr<GridPatchWrap<CoordinateSystem, Analysis> > MakeGridPatch(
         const typename CoordinateSystem::PositionVector& llc,
         const typename CoordinateSystem::PositionVector& urc,
         ConstVertexDistributionPointer vd1,
         ConstVertexDistributionPointer vd2)
      {
         return std::make_shared<
            GridPatchWrap<CoordinateSystem, Analysis> >(
               llc, urc, vd1, vd2);
      }

      template<typename CoordinateSystem>
      void DefGridPatch(const char* class_name)
      {
         using GridPatchType =
            GridPatchWrap<CoordinateSystem, AnalysisType>;
         using PositionVector = typename CoordinateSystem::PositionVector;
         using ShapePointer = typename GridPatchType::ShapePointer;
         using CurrentDensity = typename AnalysisType::CurrentDensity;
         using Current = typename AnalysisType::Current;
         using FirstCoordinate = typename CoordinateSystem::FirstCoordinate;
         using SecondCoordinate = typename CoordinateSystem::SecondCoordinate;

         void (GridPatchType::*set_load_func_1)(
            ShapePointer, CurrentDensity) = &GridPatchType::SetLoad;
         void (GridPatchType::*set_load_func_2)(
            ShapePointer, Current) = &GridPatchType::SetLoad;
         void (GridPatchType::*set_load_func_3)(
            ShapePointer, bp::object) = &GridPatchType::SetLoad;

         const EastBoundary<AnalysisType>& (GridPatchType::*eb)() const
            = &GridPatchType::east_boundary;
         const NorthBoundary<AnalysisType>& (GridPatchType::*nb)() const
            = &GridPatchType::north_boundary;
         const WestBoundary<AnalysisType>& (GridPatchType::*wb)() const
            = &GridPatchType::west_boundary;
         const SouthBoundary<AnalysisType>& (GridPatchType::*sb)() const
            = &GridPatchType::south_boundary;

#define DEF(METHOD)                                                 \
         .def(BOOST_PP_STRINGIZE(METHOD), &GridPatchType::METHOD)

         BPClassSharedPtr<GridPatchType>(
            class_name, bp::init<
            const PositionVector&, const PositionVector&,
            ConstVertexDistributionPointer, ConstVertexDistributionPointer>())
            DEF(SetMaterial)
            DEF(DeleteMaterial)
            DEF(DeleteLoad)
            DEF(Shift)
            DEF(Rotate)
            .def("SetLoad", set_load_func_3)
            .def("SetLoad", set_load_func_2)
            .def("SetLoad", set_load_func_1)
            .add_property(
               "west_boundary", bp::make_function(
                  wb, bp::return_value_policy<
                  bp::reference_existing_object>()))
            .add_property(
               "south_boundary", bp::make_function(
                  sb, bp::return_value_policy<
                  bp::reference_existing_object>()))
            .add_property(
               "east_boundary", bp::make_function(
                  eb, bp::return_value_policy<
                  bp::reference_existing_object>()))
            .add_property(
               "north_boundary", bp::make_function(
                  nb, bp::return_value_policy<
                  bp::reference_existing_object>()));

#undef DEF

         bp::implicitly_convertible<
            std::shared_ptr<GridPatchType>,
            std::shared_ptr<GridPatch<CoordinateSystem, AnalysisType> > >();
         bp::implicitly_convertible<
            std::shared_ptr<GridPatchType>,
            ConstGridPointer>();
         bp::implicitly_convertible<
            std::shared_ptr<GridPatch<CoordinateSystem, AnalysisType> >,
            ConstGridPointer>();
      }

      // --------------------
      // *** GridAssembly ***
      // --------------------

      struct GridAssemblyWrap : GridAssembly
      {
      private:
         using Base = GridAssembly;
         using GridIter = bp::stl_input_iterator<ConstGridPointer>;

      public:
         GridAssemblyWrap() : Base() {}

         GridAssemblyWrap(bp::object grids) : Base(
            GridIter(grids), GridIter())
         {
         }

         virtual ~GridAssemblyWrap()
         {
         }
      };

      // --------------------------
      // *** Boundary Condition ***
      // --------------------------

      template <typename Analysis>
      struct DirichletConditionWrap : DirichletCondition<Analysis>
      {
         using Base = DirichletCondition<Analysis>;
         using MagneticVectorPotential =
            typename Analysis::MagneticVectorPotential;
         using Time = typename Analysis::Time;
         using PositionVector = CartesianSystem::PositionVector;

         DirichletConditionWrap() : Base()
         {
         }

         explicit DirichletConditionWrap(
            MagneticVectorPotential potential)
            : Base(potential)
         {
         }

         explicit DirichletConditionWrap(bp::object potential)
         {
            bp::extract<MagneticVectorPotential>
               potential_extractor(potential);

            if (potential_extractor.check())
               Base::set_potential(potential_extractor());
            else
               Base::set_potential(
                  [potential](const PositionVector& p, Time t)
                  -> MagneticVectorPotential {
                     return bp::extract<MagneticVectorPotential>(
                        potential(p, t));
                  });
         }
      };

      // --- PathInGridGraph -----------------------------------------
      template <typename GridGraphPointer>
      struct PathInGridGraphWrap : PathInGridGraph<GridGraphPointer>
      {
         using Point = CartesianSystem::Point;
         using PointIter = bp::stl_input_iterator<Point>;
         using Base = PathInGridGraph<GridGraphPointer>;

         PathInGridGraphWrap(GridGraphPointer grid_graph,
                             bp::object points)
            : Base(grid_graph, PointIter(points), PointIter())
         {
         }

         template <typename Item, typename Iter = typename Item::Iter>
         Iter begin() const {
            return Base::template begin<Item>();
         }

         template <typename Item, typename Iter = typename Item::Iter>
         Iter end() const {
            return Base::template end<Item>();
         }
      };

      template <template <typename> class Type, typename Analysis>
      using BPClassConcreteBoundary = bp::class_<
         Type<Analysis>, bp::bases<Boundary<Analysis> >,
         boost::noncopyable>;
   }
}

BOOST_PYTHON_MODULE(MODULE_NAME)
{
   using namespace PACKAGE_NS;
   using namespace PACKAGE_NS::python;
   namespace bp = boost::python;

   // *** Boundaries ***
   using BoundaryType = Boundary<AnalysisType>;
   using BPClassBoundary =
      bp::class_<BoundaryType, boost::noncopyable>;
   BPClassBoundary("Boundary", bp::no_init)
      .def("SetBoundaryCoupling", &BoundaryType::SetBoundaryCoupling)
      .def("DeleteBoundaryCoupling", &BoundaryType::DeleteBoundaryCoupling)
      .def("SetBoundaryCondition", &BoundaryType::SetBoundaryCondition)
      .def("DeleteBoundaryCondition", &BoundaryType::DeleteBoundaryCondition)
      .def("HasBoundaryCondition", &BoundaryType::HasBoundaryCondition)
      .def("HasBoundaryCoupling", &BoundaryType::HasBoundaryCoupling);

   BPClassConcreteBoundary<EastBoundary, AnalysisType> (
      "EastBoundary",  bp::no_init);
   BPClassConcreteBoundary<NorthBoundary, AnalysisType>(
      "NorthBoundary", bp::no_init);
   BPClassConcreteBoundary<WestBoundary, AnalysisType> (
      "WestBoundary",  bp::no_init);
   BPClassConcreteBoundary<SouthBoundary, AnalysisType>(
      "SouthBoundary", bp::no_init);

   // *** Boundary conditions ***
   // ** Dirichlet conditions **
   BPClassSharedPtr<DirichletConditionWrap<AnalysisType> >(
      "DirichletCondition", bp::init<>())
      .def(bp::init<AnalysisType::MagneticVectorPotential>())
      .def(bp::init<bp::object>());
   bp::implicitly_convertible<
      std::shared_ptr<DirichletConditionWrap<AnalysisType> >,
      std::shared_ptr<DirichletCondition<AnalysisType> > >();
   bp::implicitly_convertible<
      std::shared_ptr<DirichletConditionWrap<AnalysisType> >,
      ConstBoundaryConditionPointer<AnalysisType> >();
   bp::implicitly_convertible<
      std::shared_ptr<DirichletCondition<AnalysisType> >,
      ConstBoundaryConditionPointer<AnalysisType> >();

   // ** Neumann condition **
   BPClassSharedPtr<NeumannCondition<AnalysisType> >(
      "NeumannCondition", bp::no_init)
      .def(bp::init<MagneticReluctivity, AnalysisType::MagneticFluxDensity>())
      .def(bp::init<double, AnalysisType::MagneticFluxDensity>());
   bp::implicitly_convertible<
      std::shared_ptr<NeumannCondition<AnalysisType> >,
      ConstBoundaryConditionPointer<AnalysisType> >();

   // *** Grids ***
   // ** GridAssembly **
   BPClassSharedPtr<GridAssemblyWrap>("GridAssembly", bp::no_init)
      .def(bp::init<bp::object>())
      .def("AddGrid", &GridAssemblyWrap::AddGrid);

   bp::implicitly_convertible<
      std::shared_ptr<GridAssemblyWrap >,
      std::shared_ptr<GridAssembly> >();
   bp::implicitly_convertible<
      std::shared_ptr<GridAssemblyWrap>, ConstGridPointer>();

   // ** GridPatch **
   DefGridPatch<CartesianSystem>("GridPatchCartesian");
   DefGridPatch<PolarSystem>("GridPatchPolar");

   // *** GridGraph ***
   using GridGraphType = fi::GridGraph<AnalysisType>;
   using GGPointer = fi::GridGraphPointer<AnalysisType>;
   BPClassSharedPtr<GridGraphType>(
      "GridGraph", bp::no_init)
      .add_property("potentials",
                    bp::range(&GridGraphType::begin<GridGraphType::Potentials>,
                              &GridGraphType::end<GridGraphType::Potentials>))
      .add_property("positions",
                    bp::range(&GridGraphType::begin<GridGraphType::Positions>,
                              &GridGraphType::end<GridGraphType::Positions>));

   bp::def("WriteToGraphviz", WriteToGraphviz<GGPointer>);
   bp::def("NextTimeStep", NextTimeStep<GGPointer>);

   // *** PathInGridGraph ***
   using PIGG = PathInGridGraphWrap<GGPointer>;
   BPClassSharedPtr<PIGG>(
      "PathInGridGraph", bp::init<GGPointer, bp::object>())
      .add_property("total_length", &PIGG::total_length)
      .add_property(
         "average_flux_density", &PIGG::average_flux_density)
      .add_property(
         "magnetic_voltage", &PIGG::magnetic_voltage)
      .add_property("potentials",
                    bp::range(&PIGG::begin<PIGG::Potentials>,
                              &PIGG::end<PIGG::Potentials>))
      .add_property("positions",
                    bp::range(&PIGG::begin<PIGG::Positions>,
                              &PIGG::end<PIGG::Positions>));

   bp::implicitly_convertible<
      std::shared_ptr<PIGG>,
      std::shared_ptr<PathInGridGraph<GGPointer> > >();

   // *** GridCompiler ***
   using GridCompiler = fi::GridCompiler<AnalysisType>;

   GGPointer (GridCompiler::*call_1)() = &GridCompiler::operator();
   GGPointer (GridCompiler::*call_2)(AnalysisType::Time)
      = &GridCompiler::operator();

   BPClassSharedPtr<GridCompiler>(
      "GridCompiler", bp::init<ConstGridPointer>())
      .def("__call__", call_1)
      .def("__call__", call_2);


   // Solving
   bp::def("LinearSolve", LinearSolve<GGPointer>);
   bp::def("IterativeSolve", IterativeSolve<GGPointer>, (
              bp::arg("grid_graph"),
              bp::arg("epsilon")=1e-4, bp::arg("max_iter")=36));
   bp::def("Newton", Newton<GGPointer>, (
              bp::arg("grid_graph"),
              bp::arg("epsilon")=1e-4, bp::arg("max_iter")=36));
}

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
