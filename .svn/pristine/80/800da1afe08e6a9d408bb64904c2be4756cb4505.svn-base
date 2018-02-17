###################################################################
#           Makefile for  DNS CODE - LoMa (Low Mach DNS)
###################################################################

########################################################## compiler
#FC = mpif90
FC =ifort
#FC =/usr/lib64/openmpi/bin/mpif90

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

#LFLAGS   = -O3
#FCFLAGS = -O3 -c
FCFLAGS = -c -O3 

##LFLAGS   = -O3  
##FCFLAGS = -O3 -c 

#FCFLAGS = -O3 -c -fbacktrace 
#FCFLAGS = -O3 -c -Wno-all -fbacktrace -fbounds-check
#FCFLAGS = -O3 -c -Wno-all -fbacktrace -fbounds-check
#FCFLAGS = -O0 -c
             
#CCFLAGS = -O3 -c -qmaxmem=-1 -qcache=auto -qaltivec \
#          -q64 -qarch=ppc970 -qtune=ppc970 
#CCFLAGS  = -O3 -c              
             
             
########################################################### machine
#SYSTEM = Mare Nostrum
##################################################### objects alpha
INIT = $(PWD)
SRC  = $(INIT)
OBJ  = $(INIT)/build

OBJECTS = $(OBJ)/modulos.o\
          $(OBJ)/cfdiff.o\
          $(OBJ)/rftsingle.o\
	  $(OBJ)/cftsingle.o\
          $(OBJ)/main.o 

OBJECTS2 = $(OBJ)/modulos.o\
          $(OBJ)/rftsingleclassic.o\
	  $(OBJ)/cftsingleclassic.o\
          $(OBJ)/main2.o 


############################################################ build
test: $(OBJECTS)
	@echo Linking. 
	$(FC) -o  $@ $(LFLAGS) $(OBJECTS) ../lib/libdfftpack.a 
	@echo test built, congratulations.

############################################################ compile 

############################################################ build
test2: $(OBJECTS2)
	@echo Linking. 
	$(FC) -o  $@ $(LFLAGS) $(OBJECTS2) 
	@echo test2 built, congratulations.

#################################################### compile 
#
#
$(OBJ)/modulos.o : $(SRC)/modulos.f $(SRC)/ctes3D  $(SRC)/makefile
	@echo compiling modulos.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/modulos.f

$(OBJ)/main.o : $(SRC)/main.f  $(SRC)/ctes3D $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling main.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/main.f

$(OBJ)/cfdiff.o : $(SRC)/cfdiff.f $(SRC)/ctes3D  $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling cfdiff.F
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/cfdiff.f

$(OBJ)/main2.o : $(SRC)/main2.f  $(SRC)/ctes3D $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling main2.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/main2.f

$(OBJ)/cftsingle.o : $(SRC)/cftsingle.f  $(SRC)/makefile
	@echo compiling cftsingle.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/cftsingle.f

$(OBJ)/rftsingle.o : $(SRC)/rftsingle.f  $(SRC)/makefile
	@echo compiling rftsingle.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/rftsingle.f

$(OBJ)/cftsingleclassic.o : $(SRC)/cftsingleclassic.f  $(SRC)/makefile
	@echo compiling cftsingleclassic.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/cftsingleclassic.f

$(OBJ)/rftsingleclassic.o : $(SRC)/rftsingleclassic.f  $(SRC)/makefile
	@echo compiling rftsingleclassic.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/rftsingleclassic.f


###################################################### end of file
clean:
#	find . \( -name '*.o' \) -exec rm -rf {} \;
	find $(OBJ) \( -name '*.o' \) -exec rm -rf  {} \;
	find $(OBJ) \( -name '*.mod' \) -exec rm -rf  {} \;

