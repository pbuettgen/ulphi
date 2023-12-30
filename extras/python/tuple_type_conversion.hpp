// -*- coding: utf-8 -*-
/*
 * Copyright © 2012-2016 Philipp Büttgenbach
 *
 * This file is part of ulphi, a CAE tool and library
 * for computing electromagnetic fields.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*! \file
 *
 *  \brief Conversion between C++ tuples and python tuples
 *  \author Philipp Büttgenbach
 *  \date $Date: 2016-02-14 23:53:41 +0100 (So, 14. Feb 2016) $
 *  \copyright MPL 2.0
 *
 * Conversion between C++ tuples and python tuples.  Adapted from <a
 * href="https://codereview.stackexchange.com/questions/9202/boost-python-converter-for-stdtuple">Boost
 * Python converter for std::tuple</a>.
 */

#ifndef TUPLE_CONVERSION_HPP_Wryecpen_
#define TUPLE_CONVERSION_HPP_Wryecpen_

#include <tuple>

#include <boost/python.hpp>

namespace PACKAGE_NS {
   namespace python {
      namespace detail {

         template <int ...> struct IntSequence {};

         template <int N, int ...S>
         struct GenerateIntSequence
            : GenerateIntSequence<N-1, N-1, S...> {};

         template <int ...S>
         struct GenerateIntSequence<0, S...> {
            using Type = IntSequence<S...>;
         };

         template <typename ...Args>
         struct TupleCppToPythonWrapper
         {
            using CppTuple = std::tuple<Args...>;
            using PythonTuple = boost::python::tuple;
            using IntSequenceType =
               typename GenerateIntSequence<sizeof ... (Args)>::Type;

            TupleCppToPythonWrapper(const CppTuple& in_tuple)
               : in_tuple_(in_tuple)
            {}

            PythonTuple operator()() const
            {
               return DoIt(IntSequenceType());
            }

         protected:
            template<int ...S>
            PythonTuple DoIt(IntSequence<S...>) const
            {
               return boost::python::make_tuple(std::get<S>(in_tuple_) ...);
            }

         private:
            const CppTuple& in_tuple_;
         };

         template <typename ...Args>
         struct TuplePythonToCppWrapper
         {
            using CppTuple = std::tuple<Args...>;
            using PythonTuple = boost::python::tuple;
            using IntSequenceType =
               typename GenerateIntSequence<sizeof ... (Args)>::Type;

            TuplePythonToCppWrapper(PyObject* in_tuple)
               : in_tuple_(in_tuple)
            {}

            std::tuple<Args...> operator() () const
            {
               return DoIt(IntSequenceType());
            }

         protected:
            template<int ...S>
            CppTuple DoIt(IntSequence<S...>) const
            {
               using boost::python::extract;

               return std::make_tuple(
                  (static_cast<Args>(
                     extract<Args>(PyTuple_GetItem(in_tuple_, S))))...);
            }

         private:
            PyObject* in_tuple_;
         };
      }

      template <typename ... Args>
      struct TupleToPythonConversion {
         using CppTuple = std::tuple<Args ...>;

         static PyObject* convert(const CppTuple& t)
         {
            namespace bp = boost::python;

            detail::TupleCppToPythonWrapper<Args...> wrapper(t);

            auto python_tuple = wrapper();

            return bp::incref(python_tuple.ptr());
         }
      };

      namespace detail {

         template <std::size_t i>
         struct TupleElementIsConvertible
         {
            TupleElementIsConvertible(PyObject* python_tuple)
               : python_tuple_(python_tuple)
            {}

            template <typename ... Args>
            void* operator() (std::tuple<Args...>*) const
            {
               using CppTuple = std::tuple<Args ...>;
               using boost::python::extract;
               using std::tuple_element;

               if (extract<typename tuple_element<i-1, CppTuple>::type>(
                      PyTuple_GetItem(python_tuple_, i-1)).check())
                  return TupleElementIsConvertible<i-1>(python_tuple_)(
                     (CppTuple*) nullptr);

               return nullptr;
            }

         private:
            PyObject* python_tuple_;
         };

         template <>
         struct TupleElementIsConvertible<0> {

            TupleElementIsConvertible(PyObject* python_tuple)
               : python_tuple_(python_tuple)
            {}

            template <typename ... Args>
            void* operator() (std::tuple<Args...>*)
            {
               return python_tuple_;
            }

         private:
            PyObject* python_tuple_;
         };

         template <typename ... Args>
         void* TupleCheckEachElement(
            PyObject* o, std::tuple<Args...>*)
         {
            using CppTuple = std::tuple<Args...>;

            return TupleElementIsConvertible<
               std::tuple_size<CppTuple>::value>(o)((CppTuple*) nullptr);
         }
      }

      template <typename ... Args>
      void* TupleIsConvertible(PyObject* o)
      {
         using CppTuple = std::tuple<Args ...>;

         if (not PyTuple_CheckExact(o)) return nullptr;

         Py_ssize_t tuple_size = PyTuple_Size(o);

         if (sizeof ... (Args) != tuple_size) return nullptr;

         return detail::TupleCheckEachElement(o, (CppTuple*) nullptr);
      }

      template <typename ... Args>
      void TupleConstructFromPython(
         PyObject* o,
         boost::python::converter::rvalue_from_python_stage1_data* data)
      {
         namespace bp = boost::python;

         using CppTuple = std::tuple<Args...>;

         detail::TuplePythonToCppWrapper<Args...> wrapper(o);

         void* storage =
            reinterpret_cast<
            bp::converter::rvalue_from_python_storage<CppTuple>* >(
               data)->storage.bytes;
         new (storage) CppTuple(wrapper());
         data->convertible = storage;
      }

      template <typename ... Args>
      void TupleRegisterConversion()
      {
         namespace bp = boost::python;

         using CppTuple = std::tuple<Args...>;

         // register the to-python converter
         bp::to_python_converter<
            CppTuple, TupleToPythonConversion<Args...> >();

         // register the from-python converter
         bp::converter::registry::push_back(
            &TupleIsConvertible<Args...>,
            &TupleConstructFromPython<Args...>,
            bp::type_id<CppTuple>() );
      }
   }
}

#endif

/*
 * Local Variables:
 * c-file-style: "ellemtel"
 * ispell-local-dictionary: "en_US"
 * End:
 */
