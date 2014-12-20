#  -------------------
#  Usage: 
#         gmake clean
#         gmake 
#  -------------------
#  GNU_Makefile to COMPILE and LINK  testspline.f90 and cspline.f90
#  James Orr, LSCE/IPSL, CEA-CNRS-UVSQ, 91191 Gif-sur-Yvette, France
#  15 January 2014

# To use another Fortran compiler, replace "f95" in FC and F90 lines with your compiler command
# For example, comment out 2 lines below with "f95" & uncomment the following 2 lines for "gfortran"
#FC = fort77
#FC = xlf
#FC = f95
#F90 = f95
FC = gfortran
F90 = gfortran
#FC = ifort
#F90 = ifort

#DEBUGFLAGS = -g
#LDFLAGS = -L/usr/local/lib -lnetcdf -lnetcdff
LDFLAGS = -L./ -lfastspline
INCLUDEFLAGS = -I/usr/local/include .

# List of executables to be built within the package
PROGRAMS = libfastspline.a testspline fastspline.so

# "make" builds all
all: $(PROGRAMS)

#---------------------------------------------------------------------------

SOURCES = cspline.f90 

OBJS =  cspline.o

EXEC = testspline

library = libfastspline.a
#---------------------------------------------------------------------------

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
#%: %.o
#	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General Pattern rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

#---------------------------------------------------------------------------
# Build the fastspline library containing the object files (not used, illustration only)
$(library):  $(OBJS)
	ar cr $(library) $(OBJS)

# Build the Fortran program executable that tests the fastspline library (testspline)
#$(EXEC): $(TOBJS) $(library)
#	$(FC) $(FCFLAGS) -o $(EXEC) $(TOBJS) 
$(EXEC): $(OBJS) $(library) testspline.o
	$(FC) $(FCFLAGS) -o $@ $@.f90 $(LDFLAGS)

# Build the shared object file for python
fastspline.so: $(OBJS)
	f2py -c $(SOURCES) -m fastspline --fcompiler=gnu95 --f90flags=-O3
#---------------------------------------------------------------------------

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ 

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.so

veryclean: clean
	rm -f *~ $(PROGRAMS) 
