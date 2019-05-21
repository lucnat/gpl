
UNAME_S := $(shell uname -s)

SHA1=sha1sum
FC=gfortran
FFLAGS=-fdefault-real-8 -cpp

LD=gfortran

OBJ= globals.o utils.o shuffle.o maths_functions.o mpl_module.o gpl_module.o

# Rules for main fortran files
%.o: %.f90
		@echo "F90 $@"
		@$(FC) $(FFLAGS) -c $<


# Rules for linking
test: $(OBJ) test.o
		@echo "LD $@"
		@$(LD) -o $@ $^  $(LFLAGS)

clean:
		@rm -f *.o *.mod 
