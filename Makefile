# Assumes fftw3, lapack and BLAS are installed to HOMELOCAL, 
# if not already within the paths the compiler looks for
MPIF90 = mpifort
# -Warray-temporaries
FCFLAGS = -ffpe-trap=invalid,zero,overflow -Wall -Wextra \
		  -Wno-missing-include-dirs -fimplicit-none -fexternal-blas \
		  -ffree-line-length-none -x f95-cpp-input -flto -c -O3 -m64 \
		  -I${HOMELOCAL}/include -I/usr/include 
LDFLAGS = -L${HOMELOCAL}/lib -L/usr/local \
		  -lfftw3 -lblas -llapack\
		  

# Modules
MODULES = numbers.o\
		  openmpi.o\
		  io.o\
		  parameters.o\
		  fftw.o\
		  diffops.o\
		  symmops.o\
		  test_symmops.o\
		  fieldio.o\
		  vfield.o\
		  rhs.o\
		  timestep.o\
		  stats.o\
		  lyap.o\
		  symred.o\
		  projector.o\
		  run.o\
		  solver.o

# Objects

UTILS 	= utilities/newton.o\
          utilities/eigen.o

# -------------------------------------------------------

all:  $(MODULES) main.o 
	  $(MPIF90) $(MODULES) main.o  -o dns.x $(LDFLAGS)

# -------------------------------------------------------

utils: $(MODULES) $(UTILS)
	   $(MPIF90) $(MODULES) newton.o -o newton.x $(LDFLAGS)
	   $(MPIF90) $(MODULES) eigen.o -o eigen.x $(LDFLAGS)
	   
# -------------------------------------------------------

test: $(MODULES) test.o
	  $(MPIF90) $(MODULES) test.o  -o test.x $(LDFLAGS)
	  rm -rf test 
	  mkdir test
	  cp parameters.in test
	  cp test.x test
	  cd test; ./test.x; cat d0000.txt

# -------------------------------------------------------

# compile

%.o: %.f
	$(MPIF90) $(FCFLAGS) $<

parameters.o: parameters.f90
			  bash version.sh && $(MPIF90) $(FCFLAGS) $<

%.o: %.f90
	 $(MPIF90) $(FCFLAGS) $<

clean:
	rm *.o *.mod *.x
