# Shell
SHELL		= /bin/bash

LBDR = ../../symlib/src
FOUND = false
ifeq (${F90},gfortran)  # gfortran compiler
  ifeq (${DEBUG},false)
     FFLAGS = -O3 -ffree-line-length-none -fno-underscoring -I${LBDR}
     FOUND = true
  else
     FFLAGS = -O -Wfatal-errors -fPIC -g -fbacktrace -fbounds-check -Wall -ffree-line-length-none -fno-underscoring -I${LBDR} 
     FOUND = true
  endif
endif

ifeq (${F90},ifc)  # Intel compiler
  FFLAGS = -fPIC -g -error-limit 7 -traceback -check bounds -warn  -I${LBDR}
  FOUND = true
endif

ifeq (${F90},ifort)  # Intel compiler
  ifeq (${DEBUG},false)
     FFLAGS =  -fPIC -O3 -I${LBDR} 
     FOUND = true
  else
#  F90 =  /opt/intel/fc/10.0.016/bin/ifort
     FFLAGS =  -fPIC -g -debug -error-limit 7 -heap-arrays -traceback -check bounds -warn -I${LBDR} 
     FOUND = true
#-prof-use -prof-dir .
   endif
endif


ifeq (${F90},xlf90) # IBM compiler
  FFLAGS = -g -C -fPIC -qsuffix=f=f90  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},f95) # Absoft PPC compiler
#  FFLAGS =  -profile -p ${LBDR} 
  FFLAGS = -g -Rb -Rc  -fPIC -nodefaultmod -p ${LBDR} #-ea
# B80  show entry in subprograms ; Rb bounds; Rc array conformance;
# z2 warning level
# -O3 optimization
# -ea stop after one error 
  FOUND = true
endif
ifeq (${F90},)  # If the variable isn't set, make sure compilation bombs
error:
	echo Error environmental variable "F90" not set!
endif
ifneq (${FOUND},true) # If it's set but no matching flags then bomb
error:	
	echo Error: makefile doesn\'t have flags for this compiler
endif

# if type "python" > /dev/null; then
#     python irevision.py derivative_structure_generator.f90 version.enum	
# fi

SRC = kpointgeneration.f90 wrap_kpoints.f90

OBJS = ${SRC:.f90=.o}
LIBS =  ${LBDR}/libutils.a ${LBDR}/libsym.a ${LBDR}/librational.a 

.SUFFIXES :  
.SUFFIXES :  .f .f90 .f95 .o


libkptgen.a: ${OBJS}
	ar ru $@ $?
	ranlib  $@

all: libkptgen.a kptgen.x

kptgen.x: ${OBJS} driver.o
	${F90} ${LDFLAGS} -o $@ ${OBJS} driver.o ${LIBS}

.f95.o : 
	${F90} ${FFLAGS} -c $<
.f90.o : 
	${F90} ${FFLAGS} -c $<
.f.o : 
	${F90} -c $<



CLEAN  = *.o *.mod *.a *.x 
clean : 
	rm -f ${CLEAN}
clobber : 
	rm -f  ${CLEAN}  *~ \#*
