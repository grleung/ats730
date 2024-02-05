# Makefile -- modified by Bee from Leah and Nick's example

# A simple hand-made makefile for a package including applications
# built from Fortran 90 sources, taking into account the usual
# dependency cases.

# This makefile works with the GNU make command found on
# GNU/Linux systems and often called gmake on non-GNU systems
# ======================================================================
# declarations -- modify to suit your computer, compiler, paths, etc.
# ======================================================================

# FC: The compiler
FC = pgf90
# FCFLAGS: flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -Mvect=cachesize:524288 -Munroll -Mnoframe -O2 -pc 64 -Mfree # optimization flags
FCFLAGS += -Ktrap=fp -traceback # -Mbounds # array bounds checking is slow!
# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

# libraries needed for linking, unused but kept for examples
#LDFLAGS = -li_need_this_lib

# specify paths
# MODELDIR: core model directory
MODELDIR = /home/gleung/ats730/mac_model/
# SRCDIR: directory where src code (.f90 files) are located
SRCDIR = $(MODELDIR)/src
# RUNDIR: directory from which to run the model
RUNDIR = $(MODELDIR)/run
# BINDIR: directory containing this Makefile where .o and .mod files are compiled
BINDIR = $(MODELDIR)/bin

# Executables
# PROGRAM: name of model program file (e.g. driver.f90) without the .f90 extension 
# the executable compiled in this bin directory will be named after the .f90 program file
PROGRAM = MAC
# EXE: executable name to link into run directory; can be the same as PROGRAM or not
EXE = /HW2/MAC.exe

# type "make" to compile the code (compiles all)
# after compiling, link the executable (named after the PROGRAM.f90 file) into the run 
# directory under the specified EXE name above
all: $(PROGRAM)
	ln -sf $(BINDIR)/$(PROGRAM) $(RUNDIR)/$(EXE)

# =============================================================================
# rules -- modify for your file and module .f90 names, include all dependencies
# =============================================================================

# The process of building an executable is typically done in two steps:

# 1.Compilation: every source file required for our program (typically
# x.f90 or x.F90 in Fortran) is compiled into an object file (usually x.o)
# 2.Linking: the final executable file is built by "linking" together
# all the object files compiled in the previous step; in this step,
# additional pre-compiled libraries of can be added.

model_vars.o: run_constants.o
grid.o: run_constants.o model_vars.o
mem.o: run_constants.o model_vars.o
base_state.o: thermo_functions.o grid.o constants.o run_constants.o model_vars.o
cape.o: base_state.o grid.o constants.o run_constants.o model_vars.o
io.o: model_vars.o run_constants.o 
$(PROGRAM):  io.o base_state.o grid.o constants.o run_constants.o model_vars.o cape.o thermo_functions.o mem.o

# ======================================================================
# general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files.
# Added $(SRCDIR) here in order to compile .mod and .o files in a separate 
# directory from where src code files are located
%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(SRCDIR)/%.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

# type "make clean" to remove old .mod and .o files
clean:
	rm -f *.o *.mod *.MOD

# type "make veryclean" to remove .mod and .o files as well as remove old executables
veryclean: clean
	rm -f *~ $(PROGRAM) $(RUNDIR)/$(EXE)

