###################################################################
#           Makefile for testing
###################################################################

########################################################## compiler
FC = mpif90
#FC = gfortran

#################################################### compiler flags
#COMMONBLOCKS = "cfdiff.f"
#F77FLAGS = -O3 -c -dyncom $(COMMONBLOCKS) -I/usr/lib/openmpi/include/
#F77FLAGS = -O3 -c -fast -I/usr/lib/openmpi/include/             
#LFLAGS  = -O3 -lmpi_f77 -L/usr/lib/openmpi/lib
#F77FLAGS = -O3 -c -lfftw
#LFLAGS   = -O3 
#F77FLAGS = -O3 -c
LFLAGS   = -O3 -Wl,--unresolved-symbols=ignore-in-object-files
FCFLAGS = -O3 -c 
#FCFLAGS = -O3 -c -fbacktrace 
#FCFLAGS = -O3 -c -Wno-all -fbacktrace -fbounds-check
             
#CCFLAGS = -O3 -c -qmaxmem=-1 -qcache=auto -qaltivec \
#          -q64 -qarch=ppc970 -qtune=ppc970 
#CCFLAGS  = -O3 -c              
             
             
##################################################### objects alpha
INIT = $(PWD)
SRC  = $(INIT)
OBJ  = $(INIT)
INC  =/usr/lib/openmpi/include/

OBJECTS = $(OBJ)/test.o\
	  $(OBJ)/cross.o\
	  $(OBJ)/laps.o\
	  $(OBJ)/cfdiff.o

############################################################ build
testliso: $(OBJECTS)
	@echo Linking. 
	$(FC) -o  $@ $(LFLAGS) $(OBJECTS)  #./lib/gcc/64/vmx/libfftw3f.a
	@echo liso built, congratulations.

#################################################### compile 
#
$(OBJ)/main.o : $(SRC)/test.f  $(SRC)/ctes3D $(SRC)/makefile
	@echo compiling test.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I $(INC) $(SRC) $(SRC)/test.f


$(OBJ)/cross.o : $(SRC)/cross.f  $(SRC)/ctes3D $(SRC)/makefile
	@echo compiling cross.f
	cd $(OBJ); $(FC) $(FCFLAGS)  -I$(INC) $(SRC) $(SRC)/cross.f

$(OBJ)/laps.o : $(SRC)/laps.f $(SRC)/ctes3D $(SRC)/makefile 
	@echo compiling laps.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(INC) $(SRC) $(SRC)/laps.f

$(OBJ)/cfdiff.o : $(SRC)/cfdiff.f $(SRC)/ctes3D $(SRC)/makefile
	@echo compiling cfdiff.F
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(INC) $(SRC) $(SRC)/cfdiff.f

###################################################### end of file
clean:
#	find . \( -name '*.o' \) -exec rm -rf {} \;
	find $(OBJ) \( -name '*.o' \) -exec rm -rf  {} \;

