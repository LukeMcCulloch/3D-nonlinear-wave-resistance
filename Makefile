## Simple make file, It does NOT check for interdependencies of modules.
## It probably is best to do make clean; make if you change a module.


## List main program module last
## 
SOURCES = precise.f90 \
	io.f90 \
	vtkxmlmod.f90 \
	fifthpanel.f90 \
	fourthconst.f90 \
	A.f90 \
	computeA.f90 \
	RHS.f90 \
	DynamicPanel.f90 \
	computeNonLinearFlowProperties.f90 \
	slae.f90 \
	p3.f90

## Define name of main program
PROGRAM = p3

# Compiler
FF = gfortran

# Delete program
# Linux
RM = rm -f
# DOSe
#RM = del

## Compiler options
# for the Intel Fortran 90 compiler
# CFLAGS = -c -fast -heap-arrays

# for the g95 compiler
CFLAGS = -c  -O3  -ffast-math -funroll-loops

#Linker Options
# for the Intel Fortran 90 compiler
# LDFLAGS = -fast -heap-arrays

# for the g95 compiler
#LDFLAGS = -02  -march=native -ffast-math -funroll-loops # -g # -c 03
LDFLAGS =  -march=native -ffast-math -funroll-loops # -g # -c 03


## Probably no changes necessary below this line
OBJECTS = $(SOURCES:.f90=.o)

all: $(SOURCES) $(PROGRAM)


$(PROGRAM): $(OBJECTS)
	$(FF) $(LDFLAGS) $(OBJECTS) -o $@

$(OBJECTS): %.o: %.f90
	$(FF) $(CFLAGS) $< -o $@

clean:
	$(RM) $(OBJECTS) *.mod 

realclean:
	$(RM) $(OBJECTS) *.mod $(PROGRAM)
