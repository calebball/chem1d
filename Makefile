
# I believe we use some gfortran extensions
# So it's best to be safe
FC = gfortran

# Set the compiler flags that we want to use
FASTFLAGS = -O3
DEBUGFLAGS = -g -O0 -Wall -fcheck=all -fbacktrace

# Set the paths we need
MKPATH=$(abspath $(lastword $(MAKEFILE_LIST)))
GITDIR=$(dir $(MKPATH))
SRCDIR=$(GITDIR)src
OBJDIR=$(GITDIR)build
BINDIR=$(GITDIR)bin

# Set the libraries to be used and search paths
LDFLAGS = -L/usr/local/lib -llapack -lblas -I$(OBJDIR)


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

# This makefile will be VERY unhappy if a .mod file
# is ever removed manually without its associated
# .o file also disappearing.
# If the build is not working try a `make clean`
# calling make again

chem1d: $(objfiles)
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
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/input.o : $(SRCDIR)/input.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/error.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/output.o : $(SRCDIR)/output.f90 \
	$(OBJDIR)/input.o \
	$(OBJDIR)/hartree_fock.o \
	$(OBJDIR)/moller_plesset.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/storage.o : $(SRCDIR)/storage.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/one_e_integrals.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/eri.o : $(SRCDIR)/eri.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/storage.o \
	$(OBJDIR)/true_eri_exp.o \
	$(OBJDIR)/true_eri_poly.o \
	$(OBJDIR)/quasi_eri.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/hartree_fock.o : $(SRCDIR)/hartree_fock.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/storage.o \
	$(OBJDIR)/one_e_integrals.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/moller_plesset.o : $(SRCDIR)/moller_plesset.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/storage.o \
	$(OBJDIR)/hartree_fock.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/one_e_integrals.o : $(SRCDIR)/one_e_integrals.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/special_functions.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/special_functions.o : $(SRCDIR)/special_functions.f90 \
	$(OBJDIR)/constants.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/true_eri_exp.o : $(SRCDIR)/true_eri_exp.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/special_functions.o \
	$(OBJDIR)/one_e_integrals.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/true_eri_poly.o : $(SRCDIR)/true_eri_poly.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/one_e_integrals.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/quasi_eri.o : $(SRCDIR)/quasi_eri.f90 \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/input.o \
	$(OBJDIR)/error.o \
	$(OBJDIR)/one_e_integrals.o \
	$(OBJDIR)/quasi_data.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@

$(OBJDIR)/quasi_data.o : $(SRCDIR)/quasi_data.f90 \
	$(OBJDIR)/constants.o
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@


$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@


# Rules for cleaning temporary files
.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod






