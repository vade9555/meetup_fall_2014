os := $(shell uname -s)
platform := $(shell uname -m)

ifeq ($(os), Darwin)
	F90 := gfortran
	freq := $(shell sysctl hw.cpufrequency | awk '{ print $$2 }' | cut -c 1-4)
endif
ifeq ($(os), Linux)
	F90 := ifort
	freq := $(shell cat /proc/cpuinfo | grep MHz | awk '{print $$4}' | cut -c 1-4 | uniq)
endif

experiment := $(freq)_$(os)_$(platform)

ifeq ($(F90), gfortran)
	OPTS := -Ofast -fopt-info-optimized
#	OPTS := -g -Wall -fcheck=all
endif

ifeq ($(F90), ifort)
	OPTS := -Ofast -vec-report3 -simd -xhost
endif

FFLAGS = $(OPTS)
LD = $(F90)

.SUFFIXES:
.SUFFIXES: .f90 .o

all: wrong_loop_jacobi.exe jacobi.exe blocked_jacobi.exe

blocked_jacobi.exe: blocked_jacobi.o
	$(LD) $(LDFLAGS) -o $(experiment)/$(@) $(experiment)/$^

wrong_loop_jacobi.exe: wrong_loop_jacobi.o
	$(LD) $(LDFLAGS) -o $(experiment)/$(@) $(experiment)/$^

jacobi.exe: jacobi.o
	$(LD) $(LDFLAGS) -o $(experiment)/$(@) $(experiment)/$^

.f90.o:
	-mkdir $(experiment)
	$(F90) $(FFLAGS) -c $< -o $(experiment)/$(@)

.PHONEY: clean
clean:
	$(RM) -f *.o *.exe

