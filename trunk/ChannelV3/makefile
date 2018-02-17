###################################################################
#           Makefile for  DNS CODE
###################################################################

########################################################## compiler
F77 = mpif90
#CC  = xlc_r
#F77 = ifort
#CC  = icc

#################################################### compiler flags
#COMMONBLOCKS = "cfdiff.F"
#F77FLAGS = -O3 -c -dyncom $(COMMONBLOCKS) -I/usr/lib/openmpi/include/
#F77FLAGS = -O3 -c -fast -I/usr/lib/openmpi/include/             
#LFLAGS  = -O3 -lmpi_f77 -L/usr/lib/openmpi/lib
F77FLAGS = -O3 -c -lfftw
LFLAGS   = -O3 
             
#CCFLAGS = -O3 -c -qmaxmem=-1 -qcache=auto -qaltivec \
#          -q64 -qarch=ppc970 -qtune=ppc970 
#CCFLAGS  = -O3 -c              
             
             
########################################################### machine
#SYSTEM = Mare Nostrum
##################################################### objects alpha
INIT = .
SRCDIR  = $(INIT)
OBJDIR  = $(INIT)
CALCDIR = $(INIT)

OBJECTS = $(OBJDIR)/main.o\
	  $(OBJDIR)/rftsingle.o\
	  $(OBJDIR)/cftsingle.o\
	  $(OBJDIR)/fou3D.o\
	  $(OBJDIR)/cross.o\
	  $(OBJDIR)/escru.o\
	  $(OBJDIR)/laps.o\
	  $(OBJDIR)/change.o\
	  $(OBJDIR)/cfdiff.o

############################################################ build
liso: $(OBJECTS)
	@echo Linking. 
	$(F77) -o  $@ $(LFLAGS) $(OBJECTS)  #./lib/gcc/64/vmx/libfftw3f.a
	@echo liso built, congratulations.

#################################################### compile 
#
$(OBJDIR)/main.o : $(SRCDIR)/main.F  $(SRCDIR)/ctes3D.h $(SRCDIR)/makefile
	@echo compiling main.F
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/main.F

$(OBJDIR)/cftw.o : $(SRCDIR)/cftw.F $(SRCDIR)/makefile $(SRCDIR)/ctes3D.h
	@echo compiling cftw.F
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/cftw.F

$(OBJDIR)/cftsingle.o : $(SRCDIR)/cftsingle.F $(SRCDIR)/makefile
	@echo compiling cftsingle.F
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/cftsingle.F

$(OBJDIR)/rftsingle.o : $(SRCDIR)/rftsingle.F $(SRCDIR)/makefile
	@echo compiling rftsingle.F
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/rftsingle.F

$(OBJDIR)/fou3D.o : $(SRCDIR)/fou3D.F $(SRCDIR)/ctes3D.h $(SRCDIR)/makefile
	@echo compiling fou3D.F
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/fou3D.F

$(OBJDIR)/cross.o : $(SRCDIR)/cross.F  $(SRCDIR)/ctes3D.h $(SRCDIR)/makefile
	@echo compiling cross.F
	cd $(OBJDIR); $(F77) $(F77FLAGS)  -I$(SRCDIR) $(SRCDIR)/cross.F

$(OBJDIR)/escru.o : $(SRCDIR)/escru.F  $(SRCDIR)/ctes3D.h $(SRCDIR)/makefile
	@echo compiling escru.F
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/escru.F

$(OBJDIR)/laps.o : $(SRCDIR)/laps.F $(SRCDIR)/ctes3D.h $(SRCDIR)/makefile 
	@echo compiling laps.F
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/laps.F

$(OBJDIR)/change.o : $(SRCDIR)/change.F  $(SRCDIR)/ctes3D.h $(SRCDIR)/makefile
	@echo compiling change.F
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/change.F

$(OBJDIR)/cfdiff.o : $(SRCDIR)/cfdiff.F $(SRCDIR)/ctes3D.h $(SRCDIR)/makefile
	@echo compiling cfdiff.F
	cd $(OBJDIR); $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/cfdiff.F

###################################################### end of file
clean:
	find . \( -name '*.o' \) -exec rm -rf {} \;

