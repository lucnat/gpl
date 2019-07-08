
SHELL = /bin/sh

UNAME_S := $(shell uname -s)

SHA1=sha1sum

MODE=DEBUG

FC=gfortran
AR= ar rcs
FFLAGS=-fdefault-real-8 -cpp -pedantic-errors -std=f2008
FFLAGS+= -Werror -Wall -Wno-maybe-uninitialized -Wno-uninitialized 

ifeq ($(MODE),RELEASE)
  FFLAGS += -O3 -funroll-loops -Wtaps
else
  FFLAGS += -O0 -frange-check -g -fcheck=all -Wextra
  FFLAGS += -ffpe-trap=invalid,overflow -fdump-core -fbacktrace
endif

LD=gfortran

objects=build/globals.o build/ieps.o build/utils.o build/shuffle.o build/maths_functions.o build/mpl_module.o build/gpl_module.o

libgpl.a:$(objects)
		@echo "AR $@"
		@$(AR) $@ $^

# rules to make object files into /build
build/%.o: src/%.f90
		@echo "F90 $@"
		@$(FC) $(FFLAGS) -c $< -J build -o $@

eval: libgpl.a build/eval.o
		@echo "LD $@"
		@$(LD) -o $@ build/eval.o -L. -lgpl

test: $(objects) build/test.o
		@echo "LD $@"
		@$(LD) -o $@ $^ $(LFLAGS)

clean:
		@rm -f build/*.o build/*.mod
		@rm -f test eval
