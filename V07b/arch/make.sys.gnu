# -----------------------------------------
# - CHANGE AS NECCESSARY FOR ENVIRONMENT  -
# -----------------------------------------


# -----------------------------------------
# ---- Compiler and performance flags ----
# -----------------------------------------

# GNU Compiler
FC = gfortran
FFLAGS = -O

# Intel compiler
#FC = ifort
#FFLAGS = -O2
#FFLAGS += -march=native


# -----------------------------------------
# ---- Math Libraries ----
# -----------------------------------------

# Netlib BLAS
LIBBLAS = -lblas

# Netlib LAPACK
LIBLAPACK = -llapack

# FFTW FFT
#LIBFFT = -lfftw3

# Intel OneAPI (Contains Intel MKL BLAS, LAPACK, and FFT Routines)
#LIBINTEL = -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl


# -----------------------------------------
# ---- Pre-Processor Directives ----
# -----------------------------------------
DFLAGS =


# -----------------------------------------
# ---- Parallelization Options: OpenMP ----
# -----------------------------------------
# GNU compiler
#FFLAGS += -fopenmp

# Intel compiler
#FFLAGS += -qopenmp


# -----------------------------------------
# ---- Debug Options ----
# -----------------------------------------
#FFLAGS = -g
#DFLAGS += -D__DEBUG

# Intel compiler
#FFLAGS += -warn unused,uncalled -fp-stack-check -traceback -gen-interfaces -warn interfaces -check

# Intel compiler (severe check)
#FFLAGS = -O0 -debug all -check all -warn all -extend-source 132 -traceback -gen-interfaces -fpe-all=0 -fp-stack-check -fstack-protector-all -ftrapuv -no-ftz -standard-semantics
