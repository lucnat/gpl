MODE=DEBUG
HAVE_GINAC=1


FC=gfortran
AR=ar rcs
CC=gcc
MCC=mcc
ifeq ($(HAVE_GINAC),1)
LD=g++
LFLAGS=-lgfortran
else
LD=$(FC)
endif

ifeq ($(FC),gfortran)
FFLAGS=-fdefault-real-8 -cpp -pedantic-errors -std=f2008 -J build
FFLAGS+= -Werror -Wall -Wno-maybe-uninitialized -Wno-uninitialized 

ifeq ($(MODE),RELEASE)
  FFLAGS += -O3 -march=native -mtune=native -DRELEASE
else
  FFLAGS += -O0 -frange-check -g -fcheck=all -Wextra -DDEBUG
  FFLAGS += -ffpe-trap=invalid,overflow -fdump-core -fbacktrace
endif

ifeq ($(HAVE_GINAC),1)
FFLAGS += -DHAVE_GINAC
endif

else
FFLAGS=-autodouble -module build -fpp -stand f03 -O3 -xHost -fast
endif
CFLAGS=-std=c99

files=globals.o ieps.o utils.o shuffle.o maths_functions.o mpl_module.o gpl_module.o GPL.o
objects = $(addprefix build/,$(files))

ifeq ($(FC),gfortran)
all: libgpl.a gpl eval test
else
all: libgpl.a eval test
endif

libgpl.a:$(objects)
		@echo "AR $@"
		@$(AR) $@ $^

# rules to make object files into /build
build/%.o: src/%.f90
		@echo "F90 $@"
		@$(FC) $(FFLAGS) -c $< -o $@

build/%.o: src/%.cpp
		@echo "C++ $@"
		@$(CC) -c $<  -o $@

# Mathlink related

ifeq ($(FC),gfortran)
build/mcc.internals.tmp:
		@echo "MCC --internals"
		@$(MCC) --internals > $@

build/mcc.internals:build/mcc.internals.tmp
		@echo "Base path" > $@
		@cat $@.tmp | grep "MLDK Directory" | cut -f2 -d":" >> $@
		@echo "Compiler" >> $@
		@arch=`cat $@.tmp | grep "Library Bit Type:" | cut -f2 -d":"` ; \
		 cat $@.tmp | grep "Compile Flags$$arch:" | cut -f2 -d":" >> $@
		@echo "Linker" >> $@
		@cat $@.tmp | grep "Linker Libraries" | cut -f2 -d":" >> $@

build/%.tm.c: src/%.tm build/mcc.internals
		@echo "MPREP $@"
		@$(shell sed -n '2p' build/mcc.internals)/mprep $< -o $@

build/gpl.o: build/gpl.tm.c build/mcc.internals
		@echo "CC $<"
		@$(CC) $(CFLAGS) $(shell sed -n '4p' build/mcc.internals) -o $@ -c $<

gpl: build/gpl.o libgpl.a build/mcc.internals
		@echo "LD $@"
		@$(LD)  $< libgpl.a -o $@ $(shell sed -n '6p' build/mcc.internals) -lgfortran

eval: libgpl.a build/eval.o
		@echo "LD $@"
		@$(LD) -o $@ build/eval.o -L. -lgpl $(LFLAGS)

ifeq ($(HAVE_GINAC),1)
test: $(objects) build/ginac.o build/test.o
		@echo "LD $@"
		@$(LD) -o $@ $^ -lcln -lginac $(LFLAGS)
else
test: $(objects) build/test.o
		@echo "LD $@"
		$(LD) -o $@ $^ $(LFLAGS)
endif


check: test
		./$<

clean:
		@rm -f build/*.o build/*.mod build/*.c build/mcc.internals*
		@rm -f test eval libgpl.a gpl
