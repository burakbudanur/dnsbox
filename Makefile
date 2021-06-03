# Assumes fftw3, lapack and BLAS are installed to HOMELOCAL, 
# if not already within the paths the compiler looks for
MPIF90 = mpifort
# -Warray-temporaries -ffpe-trap=invalid,zero,overflow
FCFLAGS = -Wcompare-reals -Wzerotrip -Wall -Wunused-parameter -Wdo-subscript \
          -Wno-missing-include-dirs -fimplicit-none -fexternal-blas \
		  -ffree-line-length-none -x f95-cpp-input -flto -c -O3 \
		  -I${HOMELOCAL}/include
LDFLAGS = -L${HOMELOCAL}/lib -lfftw3 -lblas -llapack

# Modules
MODULES = m_numbers.o\
		  m_openmpi.o\
		  m_io.o\
		  m_parameters.o\
		  m_work.o\
		  x_fftw.o\
		  m_fields.o\
		  m_stats.o\
		  m_state.o\
		  m_step.o\
		  m_runs.o\
		  m_solvers.o

# Objects

UTILS 	= utilities/newton.o\
          utilities/eigen.o\
		  utilities/visualize.o

# -------------------------------------------------------

all:  $(MODULES) main.o 
	  $(MPIF90) $(MODULES) main.o  -o dns.x $(LDFLAGS)

# -------------------------------------------------------

utils: $(MODULES) $(UTILS)
	   $(MPIF90) $(MODULES) newton.o -o newton.x $(LDFLAGS)
	   $(MPIF90) $(MODULES) eigen.o -o eigen.x $(LDFLAGS)
	   $(MPIF90) $(MODULES) visualize.o -o visualize.x $(LDFLAGS)
	   
# -------------------------------------------------------

# compile

%.o: %.f
	$(MPIF90) $(FCFLAGS) $<

m_parameters.o: m_parameters.f90
				bash version.sh && $(MPIF90) $(FCFLAGS) $<

%.o: %.f90
	 $(MPIF90) $(FCFLAGS) $<

clean:
	rm *.o *.mod *.x
