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
 *  \brief Grid patch class
 *  \author Philipp Büttgenbach
 *  \date $Date: 2017-06-25 00:11:26 +0200 (So, 25. Jun 2017) $
 *  \copyright MPL 2.0
 */

#ifndef GRID_PATCH_HPP_er2hvrTp_
#define GRID_PATCH_HPP_er2hvrTp_

#include <list>
#include <memory>
#include <algorithm>

#include "config.h"

#include <boost/units/io.hpp>

#include "analysis.hpp"
#include "constant_reluctivity.hpp"
#include "grid_patch_base.hpp"
#include "load_definition.hpp"
#include "material.hpp"
#include "rectangle.hpp"
#include "shape.hpp"
#include "transformation.hpp"

namespace PACKAGE_NS {
   //! A grid patch
   /*!
    * @tparam CoordinateSystem is the coordinate system type
    *
    * @tparam Analysis is the analysis type
    */
   template<typename CoordinateSystem, typename Analysis>
   struct GridPatch :
      GridPatchBase<Analysis>,
      detail::MakeShared<GridPatch<CoordinateSystem, Analysis> >
   {
   private:
      using Base = GridPatchBase<Analysis>;
      using ReturnType = typename Base::ReturnType;

   public:
      LOKI_DEFINE_CONST_VISITABLE()

      private:
      using Time = typename Analysis::Time;
      using Current = typename Analysis::Current;
      using CurrentFunction = std::function<Current (Time)>;
      using CurrentDensity = typename Analysis::CurrentDensity;
      using CurrentDensityFunction = std::function<CurrentDensity (Time)>;
      using FirstCoordinate = typename CoordinateSystem::FirstCoordinate;
      using SecondCoordinate = typename CoordinateSystem::SecondCoordinate;
      using ShapePointer = ConstShapePointer<CoordinateSystem>;
      using LoadDefPointer = ConstLoadDefinitionPointer<Analysis>;
      using LoadDescription =
         std::tuple<ShapePointer, LoadDefPointer >;
      using LoadDefinitionList = std::list<LoadDescription>;
      using MaterialDefinition =
         std::tuple<ShapePointer, ConstMaterialPointer>;
      using MaterialDefinitionList = std::list<MaterialDefinition>;
      using Point = typename CoordinateSystem::Point;
      using ShiftVector = typename CartesianSystem::Point;

   public:

      using Base::GetMaterial;
      using Base::GetLoad;

      //! Initialization
      /*!
       * @param lower_left_corner is the region's lower left corner
       *
       * @param upper_right_corner is the region's upper right corner
       *
       * @param west_east_distribution is the vertex distribution for
       * the first direction
       *
       * @param north_south_distribution is the vertex distribution for
       * the second direction
       */
      GridPatch(const Point&, const Point&, ConstVertexDistributionPointer,
                ConstVertexDistributionPointer);

      //! Initialize from an existing grid patch (unsupported)
      GridPatch(const GridPatch&) = delete;

      //! Asignment is unsupported
      GridPatch operator =(const GridPatch& other) = delete;

      //! Destruction
      virtual ~GridPatch();

      //! Assign a material to a region within this grid
      /*!
       * @param shape describes the region's shape
       * @param material is the material to assign
       */
      void SetMaterial(ShapePointer, ConstMaterialPointer);

      //! Get the material at a specified point
      /*!
       * @param point is the position at which the material shall be
       * determined.
       *
       * @return a pointer to the material if a material for the
       * specified point is found or nullpointer else.
       */
      ConstMaterialPointer GetMaterial(const Point&) const;

      ConstMaterialPointer GetMaterialHere(
         RowIndex, ColumnIndex) const override;

      ConstMaterialPointer GetMaterialHere(
         RowIndex, ColumnIndex, Axis) const override;

      //! Delete a previously made material assignment
      /*!
       * The list of material assignments is searched in reverse for a
       * shape containing the given point.
       *
       * @param point points into a material region defining shape
       */
      void DeleteMaterial(const Point&);

      //! Set a constant current density for a specified region
      /*!
       * @param shape describes the region's shape.
       *
       * @param load is the current density to apply.
       */
      void SetLoad(ShapePointer, CurrentDensity);

      //! Set a variing current density for a specified region
      /*!
       * @param shape describes the region's shape.
       *
       * @param load is a function which computes the current density
       * depending on time/frequency
       */
      void SetLoad(ShapePointer, CurrentDensityFunction);

      //! Set a constant source current for a specified region
      /*!
       * @param shape describes the region's shape.
       *
       * @param load is the total current linkage which shall be
       * applied.
       */
      void SetLoad(ShapePointer, Current);

      //! Set a time/frequency variing current
      /*!
       * @param shape describes the region's shape.
       *
       * @param load is a function which computes the current linkage
       * depending on time/frequency
       */
      void SetLoad(ShapePointer, CurrentFunction);

      //! Retrieve a previously set load condition
      /*!
       * @param point specifies for which point a load condition shall be found
       *
       * @returns the load condition.
       */
      LoadDefPointer GetLoad(const Point&) const;

      //! Delete a previously defined load
      /*!
       * All previously defined load shapes are checked in reverse
       * order whether the given point is within this shape.  As soon
       * as a shape is found, the corresponding load definition is
       * deleted.
       *
       * @param point is a point within a previously defined load
       * condition's area.
       */
      void DeleteLoad(const Point&);

      //! Shift this GridPatch
      void Shift(const ShiftVector&);

      //! Rotate this GridRegion by the given angle
      /*!
       * @param angle is the rotation angle
       */
      void Rotate(Angle);

      Length EdgeLength(RowIndex, ColumnIndex, Axis) const override;

      Length CrossingEdgeLengthHere(RowIndex, ColumnIndex, Axis) const override;

      //! Compute a vertex's position
      /*!
       * Compute a vertex's position in the grid's own (that means
       * local) coordinate system.
       *
       * @param row_index is the vertex's row index
       * @param column_index is the vertex's column index
       *
       * @return the point where the vertex is located.
       */
      Point VertexPositionLocal(RowIndex, ColumnIndex) const;

      //! Compute a point's global position
      /*!
       * Compute a point's position within the global coordinate system.
       *
       * @param point is a point within the local coordinate system.
       *
       * @returns a point within the global cartesian coordinate
       * system.
       */
      CartesianSystem::Point PositionGlobal(const Point&) const;

      //! Compute a vertex' global postition
      CartesianSystem::Point VertexPositionGlobal(
         RowIndex, ColumnIndex) const override;

      //! Compute an edge'es midpoint (local coordinates)
      /*!
       * The edge is described by its source vertex (lower left
       * corner) and its direction.
       *
       * @param row_index source vertex' row index.
       *
       * @param column_index source vertex' column index.
       *
       * @param axis direction
       */
      Point EdgeMidPointLocal(RowIndex, ColumnIndex, Axis) const;

      //! Compute an edge'es midpoint (global coordinate system)
      CartesianSystem::Point EdgeMidPointGlobal(
         RowIndex, ColumnIndex, Axis) const;

      LoadDefPointer GetLoadHere(
         RowIndex, ColumnIndex) const override;

      Area CellSizeHere(RowIndex, ColumnIndex) const override;

      Area CellSizeBottomHalf(RowIndex, ColumnIndex) const override;
      Area CellSizeTopHalf(RowIndex, ColumnIndex) const override;
      Area CellSizeLeftHalf(RowIndex, ColumnIndex) const override;
      Area CellSizeRightHalf(RowIndex, ColumnIndex) const override;

   private:
      //! Check that the corner points are meaningful
      void CheckValidCornerPoints() const;

   private:
      Point lower_left_corner_, upper_right_corner_;
      MaterialDefinitionList material_definitions_;
      LoadDefinitionList load_definitions_;
      std::list<std::unique_ptr<Transformation> > transformations_;
   };

   //--------------------------------------------------------------
   //
   //   Function implementations
   //
   //--------------------------------------------------------------

   template<typename CoordinateSystem, typename Analysis>
   GridPatch<CoordinateSystem, Analysis>::GridPatch(
      const Point& lower_left_corner, const Point& upper_right_corner,
      ConstVertexDistributionPointer west_east_distribution,
      ConstVertexDistributionPointer north_south_distribution)
      : Base(west_east_distribution, north_south_distribution),
        lower_left_corner_(lower_left_corner),
        upper_right_corner_(upper_right_corner)
   {
      using boost::units::si::constants::codata::mu_0;
      using boost::units::si::ampere_per_metre_squared;

      CheckValidCornerPoints();

      auto grid_rectangle = Rectangle<CoordinateSystem>::New(
         lower_left_corner_, upper_right_corner_);

      this->SetMaterial(
         grid_rectangle, Material::New());
      this->SetLoad(
         grid_rectangle, CurrentDensity::from_value(.0));
   }

   template<typename CoordinateSystem, typename Analysis>
   GridPatch<CoordinateSystem, Analysis>::~GridPatch()
   {
   }

   template<typename CoordinateSystem, typename Analysis>
   void GridPatch<CoordinateSystem, Analysis>::CheckValidCornerPoints() const
   {
      using std::get;

#ifdef __GNUC__
#define FUNCTION_NAME __PRETTY_FUNCTION__
#else
#define FUNCTION_NAME (typeid(*this).name())
#endif

      if (get<0>(upper_right_corner_) <= get<0>(lower_left_corner_)) {
         boost::format error_message(
            "%1%: Invalid horizontal grid dimensions %2% and %3%");
         const FirstCoordinate x0 = get<0>(lower_left_corner_);
         const FirstCoordinate x1 = get<0>(upper_right_corner_);
         throw InvalidArgument((error_message % FUNCTION_NAME % x0 % x1).str());
      }

      if (get<1>(upper_right_corner_) <= get<1>(lower_left_corner_)) {
         boost::format error_message(
            "%1%: Invalid vertical grid dimensions %2% and %3%");
         const SecondCoordinate y0 = get<1>(lower_left_corner_);
         const SecondCoordinate y1 = get<1>(upper_right_corner_);
         throw InvalidArgument((error_message % FUNCTION_NAME % y0 % y1).str());
      }
#undef FUNCTION_NAME
   }

   template<typename CoordinateSystem, typename Analysis>
   inline void GridPatch<CoordinateSystem, Analysis>::SetLoad(
      ShapePointer shape, CurrentDensity load)
   {
      load_definitions_.emplace_back(
         shape, LoadDefinition<Analysis>::New(load));
   }

   template <typename CoordinateSystem, typename Analysis>
   inline void GridPatch<CoordinateSystem, Analysis>::SetLoad(
      ShapePointer shape, CurrentDensityFunction load)
   {
      load_definitions_.emplace_back(
         shape, LoadDefinition<Analysis>::New(load));
   }

   template <typename CoordinateSystem, typename Analysis>
   inline void GridPatch<CoordinateSystem, Analysis>::SetLoad(
      ShapePointer shape, Current load)
   {
      load_definitions_.emplace_back(
         shape, LoadDefinition<Analysis>::New(load, shape->ComputeArea()));
   }

   template <typename CoordinateSystem, typename Analysis>
   inline void GridPatch<CoordinateSystem, Analysis>::SetLoad(
      ShapePointer shape, CurrentFunction load)
   {
      load_definitions_.emplace_back(
         shape, LoadDefinition<Analysis>::New(load, shape->ComputeArea()));
   }

   template <typename CoordinateSystem, typename Analysis>
   inline auto
   GridPatch<CoordinateSystem, Analysis>::GetLoadHere(
      RowIndex row_index, ColumnIndex column_index) const -> LoadDefPointer
   {
      return this->GetLoad(VertexPositionLocal(row_index, column_index));
   }

   template <typename CoordinateSystem, typename Analysis>
   inline auto
   GridPatch<CoordinateSystem, Analysis>::GetLoad(
      const Point& point) const -> LoadDefPointer
   {
      using std::get;
      using boost::units::si::ampere_per_metre_squared;

      for (auto iter = load_definitions_.rbegin(),
              end = load_definitions_.rend(); iter!=end; ++iter)
      {
         if (std::get<0>(*iter)->Within(point))
            return std::get<1>(*iter);
      }

      return LoadDefPointer(nullptr);
   }

   template<typename CoordinateSystem, typename Analysis>
   inline void GridPatch<CoordinateSystem, Analysis>::DeleteLoad(
      const Point& point)
   {
      for (auto iter = load_definitions_.rbegin(),
              end = --(load_definitions_.rend());
           iter != end; ++iter) {
         if (std::get<0>(*iter)->Within(point)) {
            load_definitions_.erase((++iter).base());
            break;
         }
      }
   }

   template<typename CoordinateSystem, typename Analysis>
   inline void GridPatch<CoordinateSystem, Analysis>::SetMaterial(
      ShapePointer shape, ConstMaterialPointer material)
   {
      material_definitions_.emplace_back(shape, material);
   }

   template<typename CoordinateSystem, typename Analysis>
   inline ConstMaterialPointer
   GridPatch<CoordinateSystem, Analysis>::GetMaterial(
      const Point& point) const
   {
      using std::get;

      for (auto material_iter = material_definitions_.rbegin(),
              material_end = material_definitions_.rend();
           material_iter != material_end; ++material_iter)
      {
         if (get<0>(*material_iter)->Within(point))
            return get<1>(*material_iter);
      }

      return ConstMaterialPointer(nullptr);
   }

   template <typename CoordinateSystem, typename Analysis>
   inline ConstMaterialPointer
   GridPatch<CoordinateSystem, Analysis>::GetMaterialHere(
      RowIndex row_index, ColumnIndex column_index) const
   {
      return this->GetMaterial(
         VertexPositionLocal(row_index, column_index));
   }

   template <typename CoordinateSystem, typename Analysis>
   inline auto GridPatch<CoordinateSystem, Analysis>::EdgeMidPointLocal(
      RowIndex row_index, ColumnIndex column_index, Axis axis) const -> Point
   {
#ifdef NDEBUG
      this->CheckValid(row_index);
      this->CheckValid(column_index);
#endif

      return ManhattenMidpoint(
         VertexPositionLocal(row_index, column_index),
         VertexPositionLocal(
            row_index + (Axis::kSecond == axis),
            column_index + (Axis::kFirst == axis)));
   }

   template <typename CoordinateSystem, typename Analysis>
   inline CartesianSystem::Point GridPatch<
      CoordinateSystem, Analysis>::EdgeMidPointGlobal(
         RowIndex row_index, ColumnIndex column_index, Axis axis) const
   {
      return this->PositionGlobal(
         this->EdgeMidPointLocal(row_index, column_index, axis));
   }

   template <typename CoordinateSystem, typename Analysis>
   inline ConstMaterialPointer
   GridPatch<CoordinateSystem, Analysis>::GetMaterialHere(
      RowIndex row_index, ColumnIndex column_index, Axis axis) const
   {
      return this->GetMaterial(
         this->EdgeMidPointLocal(row_index, column_index, axis));
   }

   template<typename CoordinateSystem, typename Analysis>
   inline void GridPatch<CoordinateSystem, Analysis>::DeleteMaterial(
      const Point& point)
   {
      for (auto iter = material_definitions_.rbegin(),
              end = --material_definitions_.rend();
           iter != end; ++iter) {
         if (std::get<0>(*iter)->Within(point)) {
            material_definitions_.erase((++iter).base());
            break;
         }
      }
   }

   template<typename CoordinateSystem, typename Analysis>
   inline auto GridPatch<CoordinateSystem, Analysis>::VertexPositionLocal(
      RowIndex row_index, ColumnIndex column_index) const -> Point
   {
      using std::get;

#ifndef NDEBUG
      this->CheckValid(row_index);
      this->CheckValid(column_index);
#endif

      return std::make_tuple(
         get<0>(lower_left_corner_)
         + (get<0>(upper_right_corner_) - get<0>(lower_left_corner_))
         * this->VertexPositionInRow(column_index),
         get<1>(lower_left_corner_)
         + (get<1>(upper_right_corner_) - get<1>(lower_left_corner_))
         * this->VertexPositionInColumn(row_index));
   }

   template<typename CoordinateSystem, typename Analysis>
   inline CartesianSystem::Point
   GridPatch<CoordinateSystem, Analysis>::PositionGlobal(
      const Point& point) const
   {
      CartesianSystem::Point cartesian_point = CoordinateSystem::ToCartesic(
         point);

      for (auto iter = transformations_.begin(), end = transformations_.end();
           iter != end; ++iter) {
         (*iter)->Apply(cartesian_point);
      }

      return cartesian_point;
   }

   template<typename CoordinateSystem, typename Analysis>
   inline CartesianSystem::Point
   GridPatch<CoordinateSystem, Analysis>::VertexPositionGlobal(
      RowIndex row_index, ColumnIndex column_index) const
   {
      return PositionGlobal(VertexPositionLocal(
                               row_index, column_index));
   }

   template<typename CoordinateSystem, typename Analysis>
   inline void GridPatch<CoordinateSystem, Analysis>::Shift(
      const ShiftVector& shift_vector)
   {
      transformations_.emplace_back(new struct Shift(shift_vector));
   }

   template<typename CoordinateSystem, typename Analysis>
   inline void GridPatch<CoordinateSystem, Analysis>::Rotate(Angle angle)
   {
      transformations_.emplace_back(new struct Rotate(angle));
   }


   template <typename CoordinateSystem, typename Analysis>
   inline Length GridPatch<
      CoordinateSystem, Analysis>::CrossingEdgeLengthHere(
         RowIndex row_index, ColumnIndex column_index, Axis axis) const
   {
#ifndef NDEBUG
      this->CheckValid(row_index);
      this->CheckValid(column_index);
      if (Axis::kFirst == axis)
         this->CheckValid(column_index + 1);
      else
         this->CheckValid(row_index + 1);
#endif

      using std::get;

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      const Point p1 = this->VertexPositionLocal(row_index, column_index);
      Point p2, p3;

      if (Axis::kFirst == axis) {
         if (RowIndex(0) == row_index) {
            p2 = ManhattenMidpoint(
               this->VertexPositionLocal(row_index+1, column_index+1), p1);
            p3 = std::make_tuple(get<0>(p2), 2.*get<1>(p1)-get<1>(p2));
         }
         else if (max_row_index == row_index) {
            p2 = ManhattenMidpoint(
               this->VertexPositionLocal(row_index-1, column_index+1), p1);
            p3 = std::make_tuple(get<0>(p2), 2.*get<1>(p1)-get<1>(p2));
         }
         else {
            p2 = ManhattenMidpoint(
               this->VertexPositionLocal(row_index+1, column_index+1), p1);
            p3 = ManhattenMidpoint(
               this->VertexPositionLocal(row_index-1, column_index+1), p1);
         }
      }
      else {
         if (ColumnIndex(0) == column_index) {
            p2 = ManhattenMidpoint(
               this->VertexPositionLocal(row_index+1, column_index+1), p1);
            p3 = std::make_tuple(2.*get<0>(p1)-get<0>(p2), get<1>(p2));
         }
         else if (max_column_index == column_index) {
            p2 = ManhattenMidpoint(
               this->VertexPositionLocal(row_index+1, column_index-1), p1);
            p3 = std::make_tuple(2.*get<0>(p1)-get<0>(p2), get<1>(p2));
         }
         else {
            p2 = ManhattenMidpoint(
               this->VertexPositionLocal(row_index+1, column_index+1), p1);
            p3 = ManhattenMidpoint(
               this->VertexPositionLocal(row_index+1, column_index-1), p1);
         }
      }

      return CoordinateSystem::ManhattenDistanceMax(p2, p3);
   }

   template<typename CoordinateSystem, typename Analysis>
   inline Length GridPatch<CoordinateSystem, Analysis>::EdgeLength(
      RowIndex row_index, ColumnIndex column_index, Axis axis) const
   {
#ifndef NDEBUG
      if (Axis::kFirst == axis)
         this->CheckValid(column_index + 1);
      else
         this->CheckValid(row_index + 1);
#endif

      return CoordinateSystem::ManhattenDistanceMax(
         this->VertexPositionLocal(row_index, column_index),
         this->VertexPositionLocal(
            row_index + (Axis::kSecond == axis),
            column_index + (Axis::kFirst == axis)));
   }

   template <typename CoordinateSystem, typename Analysis>
   inline Area GridPatch<CoordinateSystem, Analysis>::CellSizeBottomHalf(
      RowIndex row_index, ColumnIndex column_index) const
   {
      return CellSizeHere(row_index, column_index)/2.;
   }

   template <typename CoordinateSystem, typename Analysis>
   inline Area GridPatch<CoordinateSystem, Analysis>::CellSizeTopHalf(
      RowIndex row_index, ColumnIndex column_index) const
   {
      return CellSizeHere(row_index, column_index)/2.;
   }

   template <typename CoordinateSystem, typename Analysis>
   inline Area GridPatch<CoordinateSystem, Analysis>::CellSizeLeftHalf(
      RowIndex row_index, ColumnIndex column_index) const
   {
      using std::get;
      using std::make_tuple;

#ifndef NDEBUG
      this->CheckValid(row_index);
      this->CheckValid(column_index);
#endif

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      const Point this_vertex_position = VertexPositionLocal(
         row_index, column_index);

      Point lower_left_corner, upper_right_corner;

      if (0 == row_index) {
         if (0 == column_index) {
            // South west corner
            const Point p = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index+1));
            upper_right_corner = make_tuple(
               get<0>(this_vertex_position), get<1>(p));
            lower_left_corner = make_tuple(
               2.*get<0>(this_vertex_position)-get<0>(p),
               2.*get<1>(this_vertex_position)-get<1>(p));
         }
         else if (max_column_index == column_index) {
            // South east corner
            const Point upper_left_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index-1));
            lower_left_corner = make_tuple(
               get<0>(upper_left_corner),
               2.*get<1>(this_vertex_position)- get<1>(upper_left_corner));
            upper_right_corner = make_tuple(
               get<0>(this_vertex_position), get<1>(upper_left_corner));
         }
         else {
            // South boundary
            const Point upper_left_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index-1));
            upper_right_corner = make_tuple(
               get<0>(this_vertex_position), get<1>(upper_left_corner));
            lower_left_corner = make_tuple(
               get<0>(upper_left_corner),
               2.*get<1>(this_vertex_position) - get<1>(upper_left_corner));
         }
      }
      else if (max_row_index == row_index) {
         if (0 == column_index) {
            // North west corner
            const Point p = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index+1));
            lower_left_corner = make_tuple(
               2.*get<0>(this_vertex_position) - get<0>(p), get<1>(p));
            upper_right_corner = make_tuple(
               get<0>(this_vertex_position),
               2.*get<1>(this_vertex_position) - get<1>(p));
         }
         else if (max_column_index == column_index) {
            // North east corner
            lower_left_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index-1));
            upper_right_corner = make_tuple(
               get<0>(this_vertex_position),
               2.*get<1>(this_vertex_position) - get<1>(lower_left_corner));
         }
         else {
            // North boundary
            const Point lower_right_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index+1));
            lower_left_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index-1));
            upper_right_corner = make_tuple(
               std::get<0>(lower_right_corner),
               2.*get<1>(this_vertex_position)- get<1>(lower_right_corner));
         }
      }
      else if (0 == column_index) {
         // West boundary
         const Point p1 = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index-1, column_index+1));
         const Point p2 = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index+1, column_index+1));
         upper_right_corner = make_tuple(
            get<0>(this_vertex_position), get<1>(p2));
         lower_left_corner = make_tuple(
            2.*get<0>(this_vertex_position)-get<0>(p1), get<1>(p1));
      }
      else if (max_column_index == column_index) {
         // East boundary
         const Point upper_left_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index+1, column_index-1));
         lower_left_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index-1, column_index-1));
         upper_right_corner = make_tuple(
            get<0>(this_vertex_position), get<1>(upper_left_corner));
      }
      else {
         lower_left_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index-1, column_index-1));
         const Point p = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index+1, column_index+1));
         upper_right_corner = make_tuple(
            get<0>(this_vertex_position), get<1>(p));
      }

      return CoordinateSystem::SurfaceElement(
         lower_left_corner, upper_right_corner);
   }

   template <typename CoordinateSystem, typename Analysis>
   inline Area GridPatch<CoordinateSystem, Analysis>::CellSizeRightHalf(
      RowIndex row_index, ColumnIndex column_index) const
   {
      using std::get;
      using std::make_tuple;

#ifndef NDEBUG
      this->CheckValid(row_index);
      this->CheckValid(column_index);
#endif

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      const Point this_vertex_position = VertexPositionLocal(
         row_index, column_index);

      Point lower_left_corner, upper_right_corner;

      if (0 == row_index) {
         if (0 == column_index) {
            // South west corner
            upper_right_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index+1));
            lower_left_corner = make_tuple(
               get<0>(this_vertex_position),
               2.*get<1>(this_vertex_position)-get<1>(upper_right_corner));
         }
         else if (max_column_index == column_index) {
            // South east corner
            const Point p = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index-1));
            lower_left_corner = make_tuple(
               get<0>(this_vertex_position),
               2.*get<1>(this_vertex_position)- get<1>(p));
            upper_right_corner = make_tuple(
               2.*get<0>(this_vertex_position)-get<0>(p), get<1>(p));
         }
         else {
            // South boundary
            upper_right_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index+1));
            lower_left_corner = make_tuple(
               get<0>(this_vertex_position),
               2.*get<1>(this_vertex_position) - get<1>(upper_right_corner));
         }
      }
      else if (max_row_index == row_index) {
         if (0 == column_index) {
            // North west corner
            const Point lower_right_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index+1));
            lower_left_corner = make_tuple(
               get<0>(this_vertex_position), get<1>(lower_right_corner));
            upper_right_corner = make_tuple(
               get<0>(lower_right_corner),
               2.*get<1>(this_vertex_position) - get<1>(lower_right_corner));
         }
         else if (max_column_index == column_index) {
            // North east corner
            const Point p = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index-1));
            lower_left_corner = make_tuple(
               get<0>(this_vertex_position), get<1>(p));
            upper_right_corner = make_tuple(
               2.*get<0>(this_vertex_position) - get<0>(p),
               2.*get<1>(this_vertex_position) - get<1>(p));
         }
         else {
            // North boundary
            const Point lower_right_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index+1));
            lower_left_corner = make_tuple(
               get<0>(this_vertex_position), get<1>(lower_right_corner));
            upper_right_corner = make_tuple(
               std::get<0>(lower_right_corner),
               2.*get<1>(this_vertex_position)- get<1>(lower_right_corner));
         }
      }
      else if (0 == column_index) {
         // West boundary
         const Point lower_right_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index-1, column_index+1));
         upper_right_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index+1, column_index+1));
         lower_left_corner = make_tuple(
            get<0>(this_vertex_position), get<1>(lower_right_corner));
      }
      else if (max_column_index == column_index) {
         // East boundary
         const Point p2 = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index+1, column_index-1));
         const Point p1 = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index-1, column_index-1));
         lower_left_corner = make_tuple(
            get<0>(this_vertex_position), get<1>(p1));
         upper_right_corner = make_tuple(
            2.*get<0>(this_vertex_position)-get<0>(p2), get<1>(p2));
      }
      else {
         const Point p = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index-1, column_index-1));
         lower_left_corner = make_tuple(
            get<0>(this_vertex_position), get<1>(p));
         upper_right_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index+1, column_index+1));
      }

      return CoordinateSystem::SurfaceElement(
         lower_left_corner, upper_right_corner);
   }

   template <typename CoordinateSystem, typename Analysis>
   inline Area GridPatch<CoordinateSystem, Analysis>::CellSizeHere(
      RowIndex row_index, ColumnIndex column_index) const
   {
      using std::get;
      using std::make_tuple;

#ifndef NDEBUG
      this->CheckValid(row_index);
      this->CheckValid(column_index);
#endif

      const auto max_row_index = this->MaxRowIndex();
      const auto max_column_index = this->MaxColumnIndex();

      const Point this_vertex_position = VertexPositionLocal(
         row_index, column_index);

      Point lower_left_corner, upper_right_corner;

      if (0 == row_index) {
         if (0 == column_index) {
            // South west corner
            upper_right_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index+1));
            lower_left_corner = make_tuple(
               2.*get<0>(this_vertex_position)-get<0>(upper_right_corner),
               2.*get<1>(this_vertex_position)-get<1>(upper_right_corner));
         }
         else if (max_column_index == column_index) {
            // South east corner
            const Point upper_left_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index-1));
            lower_left_corner = make_tuple(
               get<0>(upper_left_corner),
               2.*get<1>(this_vertex_position)- get<1>(upper_left_corner));
            upper_right_corner = make_tuple(
               2.*get<0>(this_vertex_position)-get<0>(upper_left_corner),
               get<1>(upper_left_corner));
         }
         else {
            // South boundary
            const Point upper_left_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index-1));
            upper_right_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index+1, column_index+1));
            lower_left_corner = make_tuple(
               get<0>(upper_left_corner),
               2.*get<1>(this_vertex_position) - get<1>(upper_left_corner));
         }
      }
      else if (max_row_index == row_index) {
         if (0 == column_index) {
            // North west corner
            const Point lower_right_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index+1));
            lower_left_corner = make_tuple(
               2.*get<0>(this_vertex_position) - get<0>(lower_right_corner),
               get<1>(lower_right_corner));
            upper_right_corner = make_tuple(
               get<0>(lower_right_corner),
               2.*get<1>(this_vertex_position) - get<1>(lower_right_corner));
         }
         else if (max_column_index == column_index) {
            // North east corner
            lower_left_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index-1));
            upper_right_corner = make_tuple(
               2.*get<0>(this_vertex_position) - get<0>(lower_left_corner),
               2.*get<1>(this_vertex_position) - get<1>(lower_left_corner));
         }
         else {
            // North boundary
            const Point lower_right_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index+1));
            lower_left_corner = ManhattenMidpoint(
               this_vertex_position,
               VertexPositionLocal(row_index-1, column_index-1));
            upper_right_corner = make_tuple(
               std::get<0>(lower_right_corner),
               2.*get<1>(this_vertex_position)- get<1>(lower_right_corner));
         }
      }
      else if (0 == column_index) {
         // West boundary
         const Point lower_right_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index-1, column_index+1));
         upper_right_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index+1, column_index+1));
         lower_left_corner = make_tuple(
            2.*get<0>(this_vertex_position)-get<0>(upper_right_corner),
            get<1>(lower_right_corner));
      }
      else if (max_column_index == column_index) {
         // East boundary
         const Point upper_left_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index+1, column_index-1));
         lower_left_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index-1, column_index-1));
         upper_right_corner = make_tuple(
            2.*get<0>(this_vertex_position)-get<0>(upper_left_corner),
            get<1>(upper_left_corner));
      }
      else {
         lower_left_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index-1, column_index-1));
         upper_right_corner = ManhattenMidpoint(
            this_vertex_position,
            VertexPositionLocal(row_index+1, column_index+1));
      }

      return CoordinateSystem::SurfaceElement(
         lower_left_corner, upper_right_corner);
   }

#ifndef INSTANTIATE_
   extern template struct GridPatch<CartesianSystem, TimeTransientAnalysis<> >;
   extern template struct GridPatch<CartesianSystem, MagnetoHarmonicAnalysis<> >;
   extern template struct GridPatch<PolarSystem, TimeTransientAnalysis<> >;
   extern template struct GridPatch<PolarSystem, MagnetoHarmonicAnalysis<> >;
#endif

}

#endif // GRID_PATCH_HPP_er2hvrTp_

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
