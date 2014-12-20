fastspline
==========

Fast smoothing spline routine in Fortran 90 usable in python

A Matlab script (Weinert, 2009) was converted to Fortran 90 for
improved computational speed. The arguments are identical in both the
Fotran 90 and Matlab routines.

The associated makefile also produces a shared object file, usable in
python with the following command:

  `from fastspline import mcspline`

Then see how to use the routine with

  `print mcspline.__doc__`

REFERENCE
---------
Weinert, H. L. (2009) A fast compact algorithm for cubic spline
smoothing Computational Statistics and Data Analysis 53, 932-940.
