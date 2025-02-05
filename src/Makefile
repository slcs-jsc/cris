# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = day2doy doy2day map_pert map_rad noise_pert perturbation sampling spec2tab variance volcano

# Tests...
TESTS = pert_test

# Installation directory...
DESTDIR ?= ../bin

# Include directories...
INCDIR += -I ../libs/build/include

# Library directories...
LIBDIR += -L ../libs/build/lib 

# Static compilation...
STATIC ?= 1

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# Compile for coverage testing...
COV ?= 0

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=4 --dirty --always --tags)

# Compiler...
CC ?= gcc
MPICC ?= mpicc

# CFLAGS...
CFLAGS += $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" $(OPT) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums $(OPT) -g -fopenmp -mcmodel=medium

# LDFLAGS...
LDFLAGS += $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm
ifeq ($(STATIC),1)
  LDFLAGS +=  -lhdf5_hl -lhdf5 -lm -lz -lsz -lcurl
endif

# Optimization information...
ifeq ($(INFO),1)
  CFLAGS += -fopt-info
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -pg
endif

# Compile for coverage testing...
ifeq ($(COV),1)
  CFLAGS += --coverage
endif

# Static compilation...
ifeq ($(STATIC),1)
  CFLAGS += -static
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

.PHONY : all check clean coverage cppcheck dist doxygen indent install lizard mkdocs strip uninstall

all: $(EXC)
	rm -f *~

$(EXC): %: %.c libcris.o jurassic.o
	$(CC) $(CFLAGS) -o $@ $< libcris.o jurassic.o $(LDFLAGS)

libcris.o: libcris.c libcris.h Makefile
	$(CC) $(CFLAGS) -c -o libcris.o libcris.c

jurassic.o: jurassic.c jurassic.h Makefile
	$(CC) $(CFLAGS) -c -o jurassic.o jurassic.c

check: $(TESTS)

$(TESTS): all
	@(echo "\n===== Running \"$@\" ... =====") ; \
	  cd ../tests/$@ ; ./run.sh \
	  && (echo "\n===== Test \"$@\" passed! =====") \
	  || (echo "\n===== Test \"$@\" failed! =====" ; exit 1)

clean:
	rm -rf $(EXC) *.o *.gcda *.gcno *.gcov coverage* *~

coverage:
	lcov --capture --directory . --output-file=coverage.info ; \
	genhtml coverage.info --legend --frames --output-directory=coverage ; \
	gcov -abcfu *.c || true

cppcheck:
	cppcheck --enable=all ./

dist:
	cd ..  && tar cvfz cris_$(VERSION).tgz bin/* \
          docs/Doxyfile example/* libs/build.sh libs/*gz src/* tests/*

doxygen:
	cd ../docs && doxygen

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install:
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

lizard:
	lizard -s cyclomatic_complexity

mkdocs:
	cd ../docs && mkdocs build

strip:
	strip $(EXC)

uninstall:
	cd $(DESTDIR) && rm $(EXC)
