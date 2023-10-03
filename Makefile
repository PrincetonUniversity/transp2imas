# Makefile to link transp2imas with module IMAS library v. 3.8.1
#
# prerequisites:
#  -- a local TRANSP installation that is based on the GCC compiler

F90 = gfortran
COPTS = -cpp -g -O0 -fPIE
#COPTS = -cpp -O2

INCLUDES = -I$(IMAS_PREFIX)/include/gfortran

LIBS =	-L$(IMAS_PREFIX)/lib -limas-gfortran -limas

verify : verify.o
	@echo "   * $^ -> $@"
	$(F90) $(COPTS) $^ $(LIBS) -o $@

verify.o : verify.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

.PHONY : clean
clean :
	rm -f verify *.o *.mod *~
