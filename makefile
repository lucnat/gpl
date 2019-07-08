
SHELL = /bin/sh

UNAME_S := $(shell uname -s)

SHA1=sha1sum
FC=gfortran
FFLAGS=-fdefault-real-8 -cpp -g -pedantic-errors -Werror -std=f2008 \
       -Wall -Wno-maybe-uninitialized -Wno-uninitialized -O3 -fcheck=all

LD=gfortran

testobjects= build/globals.o build/utils.o build/shuffle.o build/maths_functions.o build/mpl_module.o build/gpl_module.o build/test.o
evaluation_objects= build/globals.o build/utils.o build/shuffle.o build/maths_functions.o build/mpl_module.o build/gpl_module.o build/eval.o

# rules to make object files into /build
build/%.o: src/%.f90
		@echo "F90 $@"
		@$(FC) $(FFLAGS) -c $< -J build -o $@

eval: $(evaluation_objects)
		@echo "LD $@"
		@$(LD) -o $@ $^ $(LFLAGS)

test: $(testobjects)
		@echo "LD $@"
		@$(LD) -o $@ $^ $(LFLAGS)

clean:
		@rm -f build/*.o build/*.mod
		@rm -f test eval
