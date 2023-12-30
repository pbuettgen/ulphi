<!-- -*- coding: utf-8 -*- -->

<!-- 
	Copyright © 2016 Philipp Büttgenbach

	All rights reserved.
-->

# Ulpi &mdash; A Python Module for Magnetostatics 


## Abstract

Ulphi is a new python module for magneto-static computations.
Contrary to most other programs of this type, it is based onto the
finite integration technique (FIT).  This choice makes it possible to
implement some special features:

* Combine cartesian and polar sections within one model
* Strong support for anisotropic materials

It can be used for all types of magneto-static and quasi-static
computations in two dimensions.


## Introduction

When having studied electrical engineering you might remember that in
the lectures it is usually postulated that the flux is distributed
homogeneously over an transformer core's cross section.  But is this
always true?  And where are the limits of this postulated
precondition?

Especially for wound core transformers there are good reasons for the
flux not being distributed homogeneously over the iron core's cross
section as the outer layer is longer than the inner layer and thereby
has a higher magnetic resistance.

To answer these questions it is necessary to model the iron core in
every detail.


## Modeling the wound iron core

A wound iron core has some special features which must be considered
within the model:

1. The magnetic resistance parallel to the windings is much lower than
   perpendicular to the windings.  This is modeled using an
   anisotropic material.
2. There are straight and rounded sections.


### Modeling with a standard FEM program

The first condition is (in some limits) fulfilled by many available
FEM programs.  But what about the second condition?  How can you model
the anisotropic material characteristic in this rounded section?  As
it turns out this is really hard because FEM programs usually support
only a cartesian coordinate system.  So, the only option is to define
the material properties section wise which is obviously only an
approximation and leads to discontinuities at the joints.  And there
is another shortcoming:  As in FEM usually a triangle (or tetrahedral)
mesh is used material properties must be interpolated for other
directions than originally given when building the model.

With all these approximations and interpolations it is really
questionable whether you get a reliable model or just a dirty hack.


## Building a clean and reliable model

In order to overcome these shortcomings a new code for magneto-statics
has been developed.  It is based onto the finite integration technique
(FIT) and supports combining cartesian and polar sections within one
model.  Thereby all shortcomings can be avoided:

* The FIT-grid strictly follows the lamination direction so that an
  interpolation of material properties becomes unnecessary.
* As polar sections are supported this is also true for the rounded
  sections where discontinuities are avoided.


## Ulphi

This new method has been implemented in a C++ library called Ulphi.
In order to avoid the time consuming compile and link cycle in
everyday use this library has been linked to python using the
[Boost python library](http://www.boost.org/doc/libs/1_61_0/libs/python/doc/html/index.html).
Postprocessing is done with [scipy](http://www.scipy.org/) and
[matplotlib](http://matplotlib.org/).


## Results

Using Ulphi only a few lines of python code are required to compute
the flux distribution over the iron core's cross section.  As it turns
out the professor was right when postulating that the flux is
homogeneously distributed over the iron core's cross section.  Only at
a really low permeability level the flux distributed becomes
inhomogeneous.


## Conclusion

The question whether the flux is distributed homogeneously over an
iron core's cross section led to the implementation of a new
magneto-static computation module called Ulphi.  It has some unique
features:

* supports cartesian and polar sections within one model
* supports anisotropic materials

It can be used for all types of magneto-static and quasi-static
computations in two dimensions.

<!--
Local Variables: 
ispell-local-dictionary: "en_US"
End:
-->

