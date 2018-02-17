###################################################################
#           Makefile for  DNS CODE - LoMa (Low Mach DNS)
###################################################################

########################################################## compiler
#F77 = mpif90
#FC =  mpixlf
FC = mpif90
#CC  = xlc_r
#F77 = ifort
#CC  = icc

#################################################### compiler flags
#COMMONBLOCKS = "cfdiff.f"
#F77FLAGS = -O3 -c -dyncom $(COMMONBLOCKS) -I/usr/lib/openmpi/include/
#F77FLAGS = -O3 -c -fast -I/usr/lib/openmpi/include/             
#LFLAGS  = -O3 -lmpi_f77 -L/usr/lib/openmpi/lib
#F77FLAGS = -O3 -c -lfftw
#LFLAGS   = -O3 
#F77FLAGS = -O3 -c
#The normal one:
#LFLAGS   = -O3 
#FCFLAGS = -c -O3 
 
########################################################### machine
#SYSTEM = MAGERIT (PWR7)

LFLAGS = -O3
FCFLAGS = -O3 -c 
#LFLAGS   = -O3 -qstrict -qcache=auto -qarch=pwr7 -qtune=pwr7
#FCFLAGS = -O3 -c -qstrict -qcache=auto -qarch=pwr7 -qtune=pwr7
#Remember to EXPORT -> ufmt_littleendian = 1-


#FCFLAGS = -O3 -c -fbacktrace 
#FCFLAGS = -O3 -c -Wno-all -fbacktrace -fbounds-check
             
#CCFLAGS = -O3 -c -qmaxmem=-1 -qcache=auto -qaltivec \
#          -q64 -qarch=ppc970 -qtune=ppc970 
#CCFLAGS  = -O3 -c              
             
             
##################################################### objects alpha
INIT = $(PWD)
SRC  = $(INIT)
OBJ  = $(INIT)/../build

OBJECTS = $(OBJ)/modulos.o\
	  $(OBJ)/rftsingle.o\
	  $(OBJ)/cftsingle.o\
	  $(OBJ)/fou3D.o\
	  $(OBJ)/cfdiff.o\
	  $(OBJ)/escru.o\
	  $(OBJ)/laps.o\
	  $(OBJ)/change.o\
	  $(OBJ)/main.o\
	  $(OBJ)/cross.o

############################################################ build
loma: $(OBJECTS)
	@echo Linking. 
	$(FC) -o  $@ $(LFLAGS) $(OBJECTS) #./lib/gcc/64/vmx/libfftw3f.a 
	@cp loma ../loma1101
	@echo loma built, congratulations.

#################################################### compile 
#

$(OBJ)/modulos.o : $(SRC)/modulos.f $(SRC)/ctes3D  $(SRC)/makefile
	@echo compiling modulos.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/modulos.f

$(OBJ)/main.o : $(SRC)/main.f  $(SRC)/ctes3D $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling main.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/main.f

$(OBJ)/cftsingle.o : $(SRC)/cftsingle.f $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling cftsingle.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/cftsingle.f

$(OBJ)/rftsingle.o : $(SRC)/rftsingle.f $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling rftsingle.F
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/rftsingle.f

$(OBJ)/fou3D.o : $(SRC)/fou3D.f $(SRC)/ctes3D $(OBJ)/modulos.o  $(SRC)/makefile
	@echo compiling fou3D.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/fou3D.f

$(OBJ)/cross.o : $(SRC)/cross.f  $(SRC)/ctes3D  $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling cross.f
	cd $(OBJ); $(FC) $(FCFLAGS)  -I$(SRC) $(SRC)/cross.f

$(OBJ)/escru.o : $(SRC)/escru.f  $(SRC)/ctes3D $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling escru.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/escru.f

$(OBJ)/laps.o : $(SRC)/laps.f $(SRC)/ctes3D $(OBJ)/modulos.o $(SRC)/makefile 
	@echo compiling laps.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/laps.f

$(OBJ)/change.o : $(SRC)/change.f  $(SRC)/ctes3D $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling change.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/change.f

$(OBJ)/cfdiff.o : $(SRC)/cfdiff.f $(SRC)/ctes3D  $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling cfdiff.F
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/cfdiff.f

###################################################### end of file
clean:
#	find . \( -name '*.o' \) -exec rm -rf {} \;
	find $(OBJ) \( -name '*.o' \) -exec rm -rf  {} \;
	find $(OBJ) \( -name '*.mod' \) -exec rm -rf  {} \;

