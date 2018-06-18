LOG		= compile.log
#FFLAGS		= -c
.SILENT:

# sources
WSRCF90		=  \
		f90wrap_generateKpoints.f90 \

SRCF90		=  \
		../src/generateKpoints.f90 \

OBJSF90		= $(SRCF90:.f90=.o)

LBD1 = ../../symlib/src
LBD2 = ../src
#Here in LIBS, we used to list needed libraries. However, the ordering
#is crucial for certain compilers and dependencies need to be listed
#in the correct order. It is more robust to just use ${LBD1}/*.a in
#the f2py-f90wrap call below; then the compiler will figure out what
#order to grab the files in for dependencies.  LIBS =
#$(LBD1)/libutils.a $(LBD1)/librational.a

#IMPORTANT: remember to include the *.o and libraries in the call to create
#the shared library! For instance, if the .o file is in ../src/*.o, then don't
#use *.o in the call, but rather ../src/*.o. It all depends on where the
#library's (i.e., the one being wrapped) .o files get compiled! If you get this
#part wrong, the .so will compile, but it will have undefined symbols for the
#subroutines you want to wrap.
all: pre_comp $(WSRCF90) $(OBJSF90) 
	-rm _kgridgen.so
	echo "Creating shared library _kgridgen.so"
	echo ${LIBS}
	-f2py-f90wrap -c -m _kgridgen $(WSRCF90) ../src/*.o ${LBD1}/*.a -I${LBD2} -I. >> /dev/null
	echo "done."

pre_comp:
	+$(MAKE) -C ../../symlib/src
	+$(MAKE) -C ../src/

$(OBJSF90): %.o: %.f90
	echo "Compiling: $^"
	-$(F90) -fPIC -I${LBD1} -I${LBD2} $^ -o $@ ${LIBS}>> $(LOG) 2>> $(LOG)
	-$(F90) -fPIC -c -I${LBD1} -I${LBD2} $^
	echo "done."

$(WSRCF90): 
	f90wrap -k kind_map -m kgridgen $(SRCF90) >> $(LOG)

clean:
	+$(MAKE) clean -C ../../symlib/src
	+$(MAKE) clean -C ../src
	-rm *.o *.mod *.so compile.log
