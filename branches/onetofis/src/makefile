###################################################################
#           Makefile for  DNS CODE
###################################################################
#
########################################################## compiler

FC= mpif90
#F77= ifort
#F77= gfortran

#################################################### compiler flags

LFLAGS = -O3
FCFLAGS = -c -O3 
#F77FLAGS = -c -O3 -tpp7 -xW
#LFLAGS = -O3 -pg -tpp7 -xW  -lm -I/usr/local/include 

##################################################### objects alpha

INIT    = $(PWD)
SRC  = $(INIT)
OBJ     = $(INIT)/../build
OBJDIR  = $(OBJ)
CALCDIR = $(INIT)

OBJECTS = $(OBJ)/modulos.o\
	  $(OBJ)/rftsingle.o\
	  $(OBJ)/cftsingle.o\
	  $(OBJ)/fou3D.o\
	  $(OBJ)/cfdiff.o\
	  $(OBJ)/escru.o\
	  $(OBJ)/laps.o\
	  $(OBJ)/change.o\
	  $(OBJ)/fou2phys.o\
 	$(OBJ)/main.o
#			 $(OBJ)/fftpack.ncar.o
	


############################################################ build

#getstats: printmsgA $(OBJECTS)
getstatsT: printmsgA $(OBJECTS)
	@echo Linking...
	$(FC) -o $@  $(OBJECTS) $(LFLAGS) 
#	$(FC) -o $@  $(OBJECTS) $(LFLAGS) /usr/lib/libfftw3f.a
	@cp getstatsT ../
	@echo getstats built, congratulations.

getstatsrho: printmsgA $(OBJECTS)
	@echo Linking...
	$(FC) -o $@  $(OBJECTS) $(LFLAGS) 
#	$(FC) -o $@  $(OBJECTS) $(LFLAGS) /usr/lib/libfftw3f.a
	@cp getstatsrho ../
	@echo getstats built, congratulations.


getstats: printmsgA $(OBJECTS)
	@echo Linking...
	$(FC) -o $@  $(OBJECTS) $(LFLAGS) 
#	$(FC) -o $@  $(OBJECTS) $(LFLAGS) /usr/lib/libfftw3f.a
	@cp getstats ../
	@echo getstats built, congratulations.

getstatsUVWscal: printmsgA $(OBJECTS)
	@echo Linking...
	$(FC) -o $@  $(OBJECTS) $(LFLAGS) 
#	$(FC) -o $@  $(OBJECTS) $(LFLAGS) /usr/lib/libfftw3f.a
	@cp getstatsUVWscal ../
	@echo getstats built, congratulations.

getstatsTscal: printmsgA $(OBJECTS)
	@echo Linking...
	$(FC) -o $@  $(OBJECTS) $(LFLAGS) 
#	$(FC) -o $@  $(OBJECTS) $(LFLAGS) /usr/lib/libfftw3f.a
	@cp getstatsTscal ../
	@echo getstats built, congratulations.






#################################################### compile 

$(OBJ)/modulos.o : $(SRC)/modulos.f $(SRC)/ctes3D  $(SRC)/makefile
	@echo compiling modulos.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/modulos.f

$(OBJ)/main.o : $(SRC)/main.f  $(SRC)/ctes3D $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling main.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/main.f

$(OBJ)/fou2phys.o : $(SRC)/main.f  $(SRC)/ctes3D $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling fou2phys.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/fou2phys.f

$(OBJ)/cftsingle.o : $(SRC)/cftsingle.f $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling cftsingle.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/cftsingle.f

$(OBJ)/rftsingle.o : $(SRC)/rftsingle.f $(OBJ)/modulos.o $(SRC)/makefile
	@echo compiling rftsingle.F
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/rftsingle.f

$(OBJ)/fou3D.o : $(SRC)/fou3D.f $(SRC)/ctes3D $(OBJ)/modulos.o  $(SRC)/makefile
	@echo compiling fou3D.f
	cd $(OBJ); $(FC) $(FCFLAGS) -I$(SRC) $(SRC)/fou3D.f

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




############################################################ build
printmsgA :
	@echo
	@echo Building for $(SYSTEM)
	@echo Compiler flags : $(FCFLAGS)
	@echo Linker   flags : $(LFLAGS)
#
####################################################### end of file
#
clean:
#	find . \( -name '*.o' \) -exec rm -rf {} \;
	find $(OBJ) \( -name '*.o' \) -exec rm -rf  {} \;
	find $(OBJ) \( -name '*.mod' \) -exec rm -rf  {} \;
