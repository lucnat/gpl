
SHELL = /bin/sh

UNAME_S := $(shell uname -s)

SHA1=sha1sum
FC=gfortran
FFLAGS=-fdefault-real-8 -cpp

LD=gfortran

testobjects= build/globals.o build/utils.o build/shuffle.o build/maths_functions.o build/mpl_module.o build/gpl_module.o build/test.o
evaluation_objects= build/globals.o build/utils.o build/shuffle.o build/maths_functions.o build/mpl_module.o build/gpl_module.o build/eval.o

# rules to make object files into /build
build/%.o: %.f90
		@echo "F90 $@"
		@$(FC) $(FFLAGS) -c $< -J build -o $@

eval: $(evaluation_objects)
		@echo "LD $@"
		@$(LD) -o $@ $^ $(LFLAGS)

test: $(testobjects)
		@echo "LD $@"
		@$(LD) -o $@ $^ $(LFLAGS)

clean:
		@rm -f *.o *.mod 
		@rm -f build/*
		@rm -f test eval
