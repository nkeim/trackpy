.. trackpy documentation master file, created by
   sphinx-quickstart on Sun Sep 16 14:53:53 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

A python implementation of Crocker-Grier
=============================================
:Release: |version|
:date: |today|

The Crocker-Grier algorithm is a method of tracking features in a
series of images from frame to frame.  The core of the algorithm is to
choose the frame-to-frame linking that globally minimizes the sum of
the squared displacements.

:mod:`trackpy` is a simple and extendable implementation of Crocker-Grier. This
version of tracky is not in pure Python; it uses
`numba <http://numba.pydata.org/>`_ to run two core algorithms at C/FORTRAN
speed, and it has gone through other design modifications to make it suitable
for tracking a large number of particles over an arbitrary number of frames.
Performance is comparable to or exceeds that of MATLAB and IDL. Note that the
`C++ implementation`__ of Crocker-Grier is still the highest-performance
option.

.. _cpp: https://github.com/tacaswell/tracking
__ cpp_


Contents:
=========

.. toctree::
   :maxdepth: 3

   reference/trackpy


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
