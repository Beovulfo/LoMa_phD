###################################################################
#           Makefile for  DNS CODE
###################################################################
#
########################################################## compiler

#F77= ifort
#F77= gfortran
FC=ifort
F77= ifort
#################################################### compiler flags

F77FLAGS = -c -O3 -vec-report
FCFLAGS = -c -O3 -vec-report
LFLAGS = -O3
#F77FLAGS = -c -O3 -tpp7 -xW
#LFLAGS = -O3 -pg -tpp7 -xW  -lm -I/usr/local/include 

##################################################### objects alpha

INIT    = $(PWD)
SRCDIR  = $(INIT)
SRC     = $(SRCDIR)
OBJ     = $(INIT)/../build
OBJDIR  = $(OBJ)
CALCDIR = $(INIT)

OBJECTS = $(OBJ)/modulos.o\
	 $(OBJ)/main.o\
	 $(OBJ)/cfdiff.o\
	 $(OBJ)/uvw.o\
	 $(OBJ)/laps.o\
	 $(OBJ)/fou3D.o\
	 $(OBJ)/swap.o\
	 $(OBJ)/escru.o\
	 $(OBJ)/cftsingle.o\
	 $(OBJ)/rftsingle.o
#			 $(OBJ)/fftpack.ncar.o
	


############################################################ build

PQRHZB: printmsgA $(OBJECTS)
	@echo Linking...
	$(F77) -o $@  $(OBJECTS) $(LFLAGS) #/usr/lib/libfftw3f.a
	@cp PQRHZB ../
	@echo tofis built, congratulations.


#################################################### compile 
$(OBJDIR)/modulos.o : $(SRCDIR)/modulos.f $(SRC)/ctes3D  $(SRC)/makefile
	@echo compiling modulos.f
	cd $(OBJ); $(F77) $(FCFLAGS) -I$(SRC) $(SRC)/modulos.f

$(OBJDIR)/main.o : $(SRCDIR)/main.f  $(SRCDIR)/ctes3D $(OBJ)/modulos.o $(SRCDIR)/makefile
	@echo compiling main.f
	cd $(OBJDIR); $(F77) $(F77FLAGS)  -I$(SRCDIR) $(SRCDIR)/main.f

$(OBJDIR)/fou3D.o : $(SRCDIR)/fou3D.f  $(SRCDIR)/ctes3D $(OBJ)/modulos.o $(SRCDIR)/makefile
	@echo compiling fou3D.f
	cd $(OBJDIR); $(F77) $(F77FLAGS)   -I$(SRCDIR) $(SRCDIR)/fou3D.f

$(OBJDIR)/cfdiff.o : $(SRCDIR)/cfdiff.f  $(SRCDIR)/ctes3D $(OBJ)/modulos.o $(SRCDIR)/makefile
	@echo compiling cfdiff.f
	cd $(OBJDIR); $(F77) $(F77FLAGS)  -I$(SRCDIR) $(SRCDIR)/cfdiff.f

$(OBJDIR)/laps.o : $(SRCDIR)/laps.f  $(SRCDIR)/ctes3D  $(OBJ)/modulos.o $(SRCDIR)/makefile
	@echo compiling laps.f
	cd $(OBJDIR); $(F77) $(F77FLAGS)    -I$(SRCDIR) $(SRCDIR)/laps.f

$(OBJDIR)/swap.o : $(SRCDIR)/swap.f  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling swap.f
	cd $(OBJDIR); $(F77) $(F77FLAGS)  -I$(SRCDIR) $(SRCDIR)/swap.f

$(OBJDIR)/escru.o : $(SRCDIR)/escru.f   $(SRCDIR)/ctes3D  $(OBJ)/modulos.o $(SRCDIR)/makefile 
	@echo compiling escru.f
	cd $(OBJDIR); $(F77) $(F77FLAGS)  -I$(SRCDIR) $(SRCDIR)/escru.f

$(OBJDIR)/uvw.o : $(SRCDIR)/uvw.f  $(SRCDIR)/ctes3D  $(OBJ)/modulos.o $(SRCDIR)/makefile
	@echo compiling uvw.f
	cd $(OBJDIR); $(F77) $(F77FLAGS)   -I$(SRCDIR) $(SRCDIR)/uvw.f

#$(OBJDIR)/cftw.o : $(SRCDIR)/cftw.f  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
#	@echo compiling cftw.f
#	cd $(OBJDIR); $(F77) $(F77FLAGS)   -I$(SRCDIR) $(SRCDIR)/cftw.f

#$(OBJDIR)/rftw.o : $(SRCDIR)/rftw.f  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
#	@echo compiling rftw.f
#	cd $(OBJDIR); $(F77) $(F77FLAGS)   -I$(SRCDIR) $(SRCDIR)/rftw.f

$(OBJDIR)/cftsingle.o : $(SRCDIR)/cftsingle.f  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling cftsingle.f
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/cftsingle.f

$(OBJDIR)/rftsingle.o : $(SRCDIR)/rftsingle.f  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
	@echo compiling rftsingle.f
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/rftsingle.f

#$(OBJDIR)/fftpack.ncar.o : $(SRCDIR)/fftpack.ncar.f  $(SRCDIR)/ctes3D $(SRCDIR)/makefile
#	@echo compiling fftpack.ncar.f
#	@cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/fftpack.ncar.f

############################################################ build
printmsgA :
	@echo
	@echo Building for $(SYSTEM)
	@echo Compiler flags : $(F77FLAGS)
	@echo Linker   flags : $(LFLAGS)
#
####################################################### end of file
#
clean:
#	find . \( -name '*.o' \) -exec rm -rf {} \;
	find $(OBJ) \( -name '*.o' \) -exec rm -rf  {} \;
	find $(OBJ) \( -name '*.mod' \) -exec rm -rf  {} \;
