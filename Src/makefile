#*****************************************************************
#*       Makefile for the MOLSIM package                         *
#*                                                               *
#*       Type: make all                                          *
#*             make ser                                          *
#*             make par                                          *
#*             make molsim_ser                                   *
#*             make molsim_par                                   *
#*             make moldyn                                       *
#*             make tidy                                         *
#*             make clean                                        *
#*             make print                                        *
#*****************************************************************

#*****************************************************************
#*                        Version                                *
#*****************************************************************

#Defines your current version from ../version.name
VER=
VERFILE = ../version.conf
ifneq ("$(wildcard $(VERFILE))","")
   VER=.$(shell cat ${VERFILE} )
   ifeq ($(strip $(shell cat ${VERFILE})),)
      VER=""
   endif
endif

BINDIR:=$(abspath ../Bin)

 
#*****************************************************************
#*                        Architecture                           *
#*****************************************************************

ARCH = LOCAL_INTEL
#ARCH = LOCAL_GFORTRAN
#ARCH = ALARIK_INTEL
#ARCH = ALARIK_GFORTRAN
#ARCH = MAC_GFORTRAN
ARCHFILE = make.arch
ifneq ("$(wildcard $(ARCHFILE))","")
   include $(ARCHFILE)
endif

#*****************************************************************
#* Use F03 c-bindings - typically dependend on library selection *
#*****************************************************************

CBIND = FALSE
CBIND = TRUE

#*****************************************************************
#*      Path to installation of FFTW library (from www.fftw.org) *
#*****************************************************************

ifeq ($(CBIND), TRUE)          #F03 c-bindings
 
  ifeq ($(ARCH), LOCAL_INTEL)
    FFTW_PATH=/home/$(USER)/.fftw
  endif
  ifeq ($(ARCH), LOCAL_GFORTRAN)
    FFTW_PATH=/home/$(USER)/.fftw
  endif
  ifeq ($(ARCH), ALARIK_INTEL)
    FFTW_PATH=/home/$(USER)/.fftw
  endif

endif

#*****************************************************************
#*                        Compile                                *
#*****************************************************************
 
#--------------
#     LOCAL_INTEL
#--------------
 
ifeq ($(ARCH), LOCAL_INTEL)
  
  COMP_SER = ifort
  LINK_SER = $(COMP_SER)


  OPT_SER  = -O3 -qopt-report -qopt-report-phase=vec   # report on vecorization
  OPT_SER  = -O0 -traceback -check bounds -check uninit -ftrapuv     # provide useful info at exceptions (use for debug)
  OPT_SER  = -O3 -ipo -fp-model strict                 # inter processoptimization
  OPT_SER  = -O3 -traceback -check all     	       # provide useful info at exceptions (gives many #406 warnings)
  OPT_SER  = -O0                                       # quick compilation (useful for debug)
  OPT_SER  = -warn     				       # provides many warnings
  OPT_SER  = -prof-gen                                 # generate profile for profiled-optimization
  OPT_SER  = -prof-use -O3 -ipo                        # profiled-optimization
  OPT_SER  = -O3                                       # normal
  OPT_SER  = -O3 -ipo                                  # inter processoptimization #proviedes approx. 15% more speed

  OPT_normal = -O3 -ipo
  OPT_quick := -O0	#quick compilation
  OPT_test := $(OPT_normal) #-fp-model strict -fp-model source
  OPT_debug := -O0 -traceback -check bounds -check uninit -ftrapuv     # provide useful info at exceptions (use for debug)
  OPT_warn := -warn     				       # provides many warnings

  mode = default
  ifeq ($(mode),normal)
     	OPT_SER = $(OPT_normal)
  endif
  ifeq ($(mode),test)
     	OPT_SER = $(OPT_test)
  endif
  ifeq ($(mode),debug)
     	OPT_SER = $(OPT_debug)
  endif
  ifeq ($(mode),quick)
     	OPT_SER = $(OPT_quick)
  endif
  ifeq ($(mode),warn)
     	OPT_SER = $(OPT_warn)
  endif
  
  COMP_PAR = $(MPIFC)

  LINK_PAR = $(COMP_PAR)
  OPT_PAR  = $(OPT_SER)
   
   LLIB =
  
endif

#--------------
#    LOCAL_GFORTRAN
#--------------
 
ifeq ($(ARCH), LOCAL_GFORTRAN)
  
  COMP_SER = gfortran
  LINK_SER = $(COMP_SER)
  OPT_SER  = -O3 -ffree-form -ffree-line-length-none
   
  LLIB =
  
endif

#--------------
#     ALARIK_INTEL
#--------------

ifeq ($(ARCH), ALARIK_INTEL)

  COMP_SER = ifort
  LINK_SER = $(COMP_SER)
  OPT_SER  = -O3

  COMP_PAR = $(MPIFC)
  LINK_PAR = $(COMP_PAR)
  OPT_PAR  = $(OPT_SER)
 
  LLIB =

# for Intel compiler on Alarik we need static linking for FFTW3 to work
#    LINK_SER += -Bstatic
#    LINK_PAR += -Bstatic

endif  

#--------------
#     ALARIK_GFORTRAN  (Keep ARCH = ALARIK_INTEL for NOW)
#--------------

ifeq ($(ARCH), ALARIK_GFORTRAN)

  COMP_SER = gfortran
  LINK_SER = $(COMP_SER)
  OPT_SER  = -O3 -ffree-form -ffree-line-length-none  -g

  COMP_PAR = $(MPIFC)
  LINK_PAR = $(COMP_PAR)
  OPT_PAR  = $(OPT_SER)
 
  LLIB =

endif

#--------------
#     MAC_GFORTRAN
#--------------
 
ifeq ($(ARCH), MAC_GFORTRAN)
  
  COMP_SER = gfortran
  LINK_SER = $(COMP_SER)
  OPT_SER  = -O3 -ffree-form -ffree-line-length-none
   
  LLIB =
  
endif
#-----------------------------------------------------------------
#                          Defines
#-----------------------------------------------------------------

DEFINES_BASE = -D$(ARCH)

ifeq ($(CBIND), TRUE)                      #F03 c-bindings
  DEFINES_BASE += -DF03_CBIND
    ifeq ($(mode),test)
      DEFINES_BASE += -D_TEST_
    endif
endif


DEFINES_SER  = $(DEFINES_BASE)
DEFINES_PAR  = $(DEFINES_BASE) -D_PAR_
 
#-----------------------------------------------------------------
#                          Compiler flags
#-----------------------------------------------------------------
 
FLAG_SER  = $(DEFINES_SER) $(OPT_SER)  
FLAG_PAR  = $(DEFINES_PAR) $(OPT_PAR)

ifeq ($(ARCH), LOCAL_INTEL)
 FLAG_SER += -I$(FFTW_PATH)/include
 FLAG_PAR += -I$(FFTW_PATH)/include
endif
ifeq ($(ARCH), LOCAL_GFORTRAN)
 FLAG_SER += -I$(FFTW_PATH)/include
 FLAG_PAR += -I$(FFTW_PATH)/include
endif
ifeq ($(ARCH), ALARIK_INTEL)
 FLAG_SER += -I$(FFTW3_HOME)/include
 FLAG_PAR += -I$(FFTW3_HOME)/include
endif
ifeq ($(ARCH), ALARIK_GFORTRAN)
 FLAG_SER += -I$(FFTW3_HOME)/include
 FLAG_PAR += -I$(FFTW3_HOME)/include
endif

#-----------------------------------------------------------------
#                          Common libraries
#-----------------------------------------------------------------


ifeq ($(ARCH), LOCAL_INTEL)
   LLIB += -L$(FFTW_PATH)/lib -lfftw3 -lm
endif
ifeq ($(ARCH), LOCAL_GFORTRAN)
   LLIB += -L$(FFTW_PATH)/lib -lfftw3 -lm
endif
ifeq ($(ARCH), ALARIK_INTEL)
   LLIB += $(FFTW3_HOME)/lib/libfftw3.so  #linker issues with ifort on Alarik => workaround
endif


#-----------------------------------------------------------------
 
OBJS1=mol.o           \
      particle.o      \
      potential.o     \
      coordinate.o    \
      md.o            \
      mc.o            \
      bd.o            \
      nlist.o         \
      energy.o        \
      denergy.o       \
      dump.o          \
      group.o         \
      static.o        \
      dynamic.o       \
      image.o         \
      statistics.o    \
      mixed.o         \
      molaux.o        \
      mollib.o        \
      moluser.o       \
      sso.o	      \
      mesh.o
 
OBJS2=mol.o           \
      particle_par.o  \
      potential_par.o \
      coordinate_par.o\
      md_par.o        \
      mc_par.o        \
      bd_par.o        \
      nlist_par.o     \
      energy_par.o    \
      denergy_par.o   \
      dump_par.o      \
      group_par.o     \
      static_par.o    \
      dynamic_par.o   \
      image.o         \
      statistics_par.o\
      mixed_par.o     \
      molaux.o        \
      mollib_par.o    \
      parallel_par.o  \
      moluser_par.o   \
      sso_par.o       \
      mesh_par.o
 
OBJS3=mol.o           \
      statistics.o    \
      molaux.o        \
      mollib.o        \
      mesh.o
 
default : def
def     :
	@echo Specify any one of: all, ser, par, molsim_ser, molsim_par, moldyn
all     : molsim_ser molsim_par
ser     : molsim_ser
par     : molsim_par
ver	:
	@echo Version is specified as \"$(VER)\"
	@echo bindir is \"$(BINDIR)\"
 
molsim_ser     : molsim.o $(OBJS1)
	$(LINK_SER) $(FLAG_SER) -o molsim_ser.exe molsim.o $(OBJS1) $(LLIB)
	mkdir -p ../Bin; mv molsim_ser.exe ../Bin/molsim_ser.exe; 
	cp molsim_ser $(HOME)/bin/molsim_ser$(VER); sed -i "s|^bin=.*|bin=$(BINDIR)|g" $(HOME)/bin/molsim_ser$(VER)
molsim_par     : molsim_par.o $(OBJS2)
	$(LINK_PAR) $(FLAG_PAR) -o molsim_par.exe molsim_par.o $(OBJS2) $(LLIB)
	mkdir -p ../Bin; mv molsim_par.exe ../Bin/molsim_par.exe; 
	cp molsim_par $(HOME)/bin/molsim_par$(VER); sed -i "s|^bin=.*|bin=$(BINDIR)|g" $(HOME)/bin/molsim_par$(VER)
moldyn         : moldyn.o $(OBJS3)
	$(COMP_SER) $(FLAG_SER) -o moldyn.exe moldyn.o $(OBJS3) $(LLIB)
	mkdir -p ../Bin; mv moldyn.exe ../Bin/moldyn.exe; cp moldyn $(HOME)/bin/moldyn$(VER)

testmesh       : testmesh.o mesh.o
	$(COMP_SER) $(FLAG_SER) -o testmesh.exe testmesh.o mesh.o
	./testmesh.exe

mol.o           : mol.F90         statistics.o mesh.o
	$(COMP_SER) $(FLAG_SER) -c mol.F90
molsim.o        : molsim.F90      mol.o      statistics.o
	$(COMP_SER) $(FLAG_SER) -c molsim.F90
moldyn.o        : moldyn.F90	statistics.o
	$(COMP_SER) $(FLAG_SER) -c moldyn.F90
particle.o      : particle.F90    mol.o mollib.o
	$(COMP_SER) $(FLAG_SER) -c particle.F90
potential.o     : potential.F90   mol.o
	$(COMP_SER) $(FLAG_SER) -c potential.F90
coordinate.o    : coordinate.F90  mol.o
	$(COMP_SER) $(FLAG_SER) -c coordinate.F90
md.o            : md.F90        mol.o
	$(COMP_SER) $(FLAG_SER) -c md.F90
mc.o            : mc.F90        mol.o	nlist.o
	$(COMP_SER) $(FLAG_SER) -c mc.F90
bd.o            : bd.F90      mol.o	md.o
	$(COMP_SER) $(FLAG_SER) -c bd.F90
nlist.o         : nlist.F90       mol.o
	$(COMP_SER) $(FLAG_SER) -c nlist.F90
energy.o        : energy.F90      mol.o	mollib.o
	$(COMP_SER) $(FLAG_SER) -c energy.F90
denergy.o       : denergy.F90     mol.o	energy.o	mollib.o
	$(COMP_SER) $(FLAG_SER) -c denergy.F90
dump.o		: dump.F90 	mol.o
	$(COMP_SER) $(FLAG_SER) -c dump.F90
group.o		: group.F90	mol.o
	$(COMP_SER) $(FLAG_SER) -c group.F90
static.o        : static.F90      mol.o      statistics.o
	$(COMP_SER) $(FLAG_SER) -c static.F90
dynamic.o        : dynamic.F90      mol.o      statistics.o
	$(COMP_SER) $(FLAG_SER) -c dynamic.F90
image.o         : image.F90      mol.o
	$(COMP_SER) $(FLAG_SER) -c image.F90
statistics.o    : statistics.F90
	$(COMP_SER) $(FLAG_SER) -c statistics.F90
mixed.o        : mixed.F90      mol.o
	$(COMP_SER) $(FLAG_SER) -c mixed.F90
molaux.o        : molaux.F90      mol.o
	$(COMP_SER) $(FLAG_SER) -c molaux.F90
mollib.o        : mollib.F90
	$(COMP_SER) $(FLAG_SER) -c mollib.F90
moluser.o       : moluser.F90  mol.o   	   statistics.o	potential.o
	$(COMP_SER) $(FLAG_SER) -c moluser.F90
mesh.o          : mesh.F90
	$(COMP_SER) $(FLAG_SER) -c mesh.F90
sso.o      : sso.F90  mc.o	mol.o
	$(COMP_SER) $(FLAG_SER) -c sso.F90
testmesh.o      : testmesh.F90  mesh.o
	$(COMP_SER) $(FLAG_SER) -c testmesh.F90
molsim_par.o    : molsim.F90      mol.o      statistics_par.o
	$(COMP_PAR) $(FLAG_PAR) -o molsim_par.o -c molsim.F90
particle_par.o  : particle.F90    mol.o mollib.o
	$(COMP_PAR) $(FLAG_PAR) -o particle_par.o -c particle.F90
potential_par.o : potential.F90   mol.o
	$(COMP_PAR) $(FLAG_PAR) -o potential_par.o -c potential.F90
coordinate_par.o : coordinate.F90 mol.o
	$(COMP_PAR) $(FLAG_PAR) -o coordinate_par.o -c coordinate.F90
md_par.o        : md.F90      mol.o
	$(COMP_PAR) $(FLAG_PAR) -o md_par.o -c md.F90
mc_par.o        : mc.F90      mol.o	nlist.o
	$(COMP_PAR) $(FLAG_PAR) -o mc_par.o -c mc.F90
bd_par.o        : bd.F90      mol.o	md_par.o
	$(COMP_PAR) $(FLAG_PAR) -o bd_par.o -c bd.F90
nlist_par.o     : nlist.F90       mol.o
	$(COMP_PAR) $(FLAG_PAR) -o nlist_par.o -c nlist.F90
energy_par.o    : energy.F90      mol.o	mollib_par.o
	$(COMP_PAR) $(FLAG_PAR) -o energy_par.o -c energy.F90
denergy_par.o   : denergy.F90 	mol.o	energy_par.o	mollib_par.o
	$(COMP_PAR) $(FLAG_PAR) -o denergy_par.o -c denergy.F90
dump_par.o	: dump.F90 	mol.o
	$(COMP_PAR) $(FLAG_PAR) -o dump_par.o -c dump.F90
group_par.o	: group.F90	mol.o
	$(COMP_PAR) $(FLAG_PAR) -o group_par.o -c group.F90
static_par.o    : static.F90      mol.o      statistics_par.o
	$(COMP_PAR) $(FLAG_PAR) -o static_par.o -c static.F90
dynamic_par.o    : dynamic.F90      mol.o      statistics_par.o
	$(COMP_PAR) $(FLAG_PAR) -o dynamic_par.o -c dynamic.F90
statistics_par.o   : statistics.F90
	$(COMP_PAR) $(FLAG_PAR) -o statistics_par.o -c statistics.F90
mixed_par.o     : mixed.F90      mol.o
	$(COMP_PAR) $(FLAG_PAR) -o mixed_par.o -c mixed.F90
mollib_par.o    : mollib.F90      mol.o
	$(COMP_PAR) $(FLAG_PAR) -o mollib_par.o -c mollib.F90
parallel_par.o  : parallel.F90    mol.o
	$(COMP_PAR) $(FLAG_PAR) -o parallel_par.o -c parallel.F90
moluser_par.o : moluser.F90  mol.o     potential_par.o    statistics_par.o
	$(COMP_PAR) $(FLAG_PAR) -o moluser_par.o -c moluser.F90
mesh_par.o    : mesh.F90
	$(COMP_PAR) $(FLAG_PAR) -o mesh_par.o -c mesh.F90
sso_par.o      : sso.F90  mc_par.o	mol.o
	$(COMP_PAR) $(FLAG_PAR) -o sso_par.o -c sso.F90
 
#*****************************************************************
#*                         Tidy                                  *
#*****************************************************************

tidy 	:
	rm -f *.o *.mod
 
#*****************************************************************
#*                         Clean                                 *
#*****************************************************************

clean 	:
	rm -f *.o *.mod *.i90 ../bin/*.exe *__genmod*
 
#*****************************************************************
#*                         Print                                 *
#*****************************************************************
 
PrintCom = lw
print	:
	$(PrintCom) *90
	$(PrintCom) molsim.ver
	$(PrintCom) molsim.lib
	$(PrintCom) molsim_ser
	$(PrintCom) molsim_par
	$(PrintCom) moldyn
	$(PrintCom) makefile
	$(PrintCom) archive
