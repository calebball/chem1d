
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

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c $^ -o $@


# Rules for cleaning temporary files
.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod






