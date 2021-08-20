# -----------------------------------------------------------
# MAKEFILE FOR SINGLE INCLUSIVE HADRON PRODUCTION (UNPOLARIZED)
# -----------------------------------------------------------
#
#
# SELECT FORTRAN COMPILER AND COMPILER OPTIONS
# 
FF = gfortran
FFLAGS= -O3
#
# OPTION 2: INTEL COMPILER (enable /psi_devel/intel-5.0 with en2)
#           (use ifc as COMPILER)
#
#FORT = ifc
#FFLAGS= -O3 -w
#
#
# COMPILE AND LINK (EXECUTABLE FILE: hadunp.x)
#
all: hadrive-ms.o \
	intvegas.o \
	grids.o \
	part-unp-ms.o \
	me-unp-ms-new.o \
	pkhff.o \
	kkp.o \
	Ctq6Pdf.o \
	Ctq5Par.o \
	alphas.o
	$(FF) $(FFLAGS) \
	hadrive-ms.o \
	intvegas.o \
	grids.o \
	part-unp-ms.o \
	me-unp-ms-new.o \
	pkhff.o \
	kkp.o \
	Ctq6Pdf.o \
	Ctq5Par.o \
	alphas.o \
	-o hadunp.x \
	$(shell lhapdf-config --ldflags)

%.o: %.f
	            $(FF) $(FFLAGS) -c $<


alphas.o: NNPDF_driver/alphas.f
	        $(FF) $(FFLAGS) -c $<

#
# CLEAN-UP
#

clean:
	rm -f *.o NNPDF_driver/*.o core *~ *.x
