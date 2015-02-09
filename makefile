
# Check for a compiler choice, otherwise use gfortran
FC = gfortran

# Set the compiler flags that we want to use
FASTFLAGS = -O3
DEBUGFLAGS = -g -O0 -Wall -fcheck=all -fbacktrace

# Set the paths we need
MKPATH=$(abspath $(lastword $(MAKEFILE_LIST)))
BASEDIR=$(dir $(MKPATH))
SRCDIR=$(BASEDIR)src
OBJDIR=$(BASEDIR)build
BINDIR=$(BASEDIR)bin

ifeq ($(strip $(INSTALL_PATH)),)
INSTALL_PATH=/usr/local/bin
endif

# Set the libraries to be used and search paths
ifeq ($(strip $(BLAS_PATH)),)
BLAS_PATH = /usr/local/lib
endif

ifeq ($(strip $(LAPACK_PATH)),)
LAPACK_PATH = /usr/local/lib
endif

LDFLAGS += -L$(BLAS_PATH) -L$(LAPACK_PATH) \
	-llapack -lblas -I$(OBJDIR)

# Targets
all: FCFLAGS=$(FASTFLAGS)
all: chem1d

debug: FCFLAGS=$(DEBUGFLAGS)
debug: chem1d


# Build all the source files in SRCDIR
srcfiles=$(wildcard $(SRCDIR)/*.f90)
modfiles=$(addprefix $(OBJDIR)/, $(addsuffix .mod, \
	$(basename $(notdir $(srcfiles)))))
objfiles=$(addprefix $(OBJDIR)/, $(addsuffix .o, \
	$(basename $(notdir $(srcfiles)))))


# Build rules
chem1d: $(objfiles)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FCFLAGS) -o $(BINDIR)/chem1d \
	$(objfiles) $(LDFLAGS)


$(OBJDIR)/chem1d.o : $(SRCDIR)/chem1d.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/output.o \
	$(OBJDIR)/storage.o \
	$(OBJDIR)/eri.o \
	$(OBJDIR)/hartree_fock.o \
	$(OBJDIR)/moller_plesset.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/input.o : $(SRCDIR)/input.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/error.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/output.o : $(SRCDIR)/output.f90 \
	$(OBJDIR)/input.o \
	$(OBJDIR)/hartree_fock.o \
	$(OBJDIR)/moller_plesset.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/storage.o : $(SRCDIR)/storage.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/one_e_integrals.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/eri.o : $(SRCDIR)/eri.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/storage.o \
	$(OBJDIR)/true_eri_exp.o \
	$(OBJDIR)/true_eri_poly.o \
	$(OBJDIR)/quasi_eri.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/hartree_fock.o : $(SRCDIR)/hartree_fock.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/storage.o \
	$(OBJDIR)/one_e_integrals.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/moller_plesset.o : $(SRCDIR)/moller_plesset.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/storage.o \
	$(OBJDIR)/hartree_fock.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/one_e_integrals.o : $(SRCDIR)/one_e_integrals.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/special_functions.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/special_functions.o : $(SRCDIR)/special_functions.f90 \
	$(OBJDIR)/constants.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/true_eri_exp.o : $(SRCDIR)/true_eri_exp.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/special_functions.o \
	$(OBJDIR)/one_e_integrals.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/true_eri_poly.o : $(SRCDIR)/true_eri_poly.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/one_e_integrals.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/quasi_eri.o : $(SRCDIR)/quasi_eri.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/one_e_integrals.o \
	$(OBJDIR)/quasi_data.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/quasi_data.o : $(SRCDIR)/quasi_data.f90 \
	$(OBJDIR)/constants.o
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $< -o $@


$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@


# Rules for cleaning build directories and installing
# the software
.PHONY: clean install

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod

install:
	cp $(BINDIR)/chem1d $(INSTALL_PATH)

