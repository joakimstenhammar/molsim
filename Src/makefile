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
    FFTW_PATH=/usr/local/include
  endif
  ifeq ($(ARCH), LOCAL_GFORTRAN)
    FFTW_PATH=/usr/local/include
  endif
  ifeq ($(ARCH), ALARIK_INTEL)
    FFTW_PATH=/$(FFTW3_HOME)/include
  endif
  ifeq ($(ARCH), ALARIK_GFORTRAN)
    FFTW_PATH=/$(FFTW3_HOME)/include
  endif

  FFTWLIB = none
  FFTWFILE = make.fftwpath
  ifneq ("$(wildcard $(FFTWFILE))","")
    include $(FFTWFILE)
  endif

endif

#*****************************************************************
#*                        Compile                                *
#*****************************************************************

#----------------
#     LOCAL_INTEL
#----------------

ifeq ($(ARCH), LOCAL_INTEL)

  COMP_SER = ifort
  LINK_SER = $(COMP_SER)


  OPT_SER = -O3 -qopt-report -qopt-report-phase=vec             # report on vecorization
  OPT_SER = -O0 -traceback -check bounds -check uninit -ftrapuv # provide useful info at exceptions (use for debug)
  OPT_SER = -O3 -ipo -fp-model strict                           # inter processoptimization
  OPT_SER = -O3 -traceback -check all                           # provide useful info at exceptions (gives many #406 warnings)
  OPT_SER = -O0                                                 # quick compilation (useful for debug)
  OPT_SER = -warn                                               # provides many warnings
  OPT_SER = -prof-gen                                           # generate profile for profiled-optimization
  OPT_SER = -prof-use -O3 -ipo                                  # profiled-optimization
  OPT_SER = -O3                                                 # normal
  OPT_SER = -O3 -ipo                                            # inter processoptimization #proviedes approx. 15% more speed

  # default compilation
  OPT_normal := -O3 -ipo
  # quick compilation:
  OPT_quick  := -O0
  # -fp-model strict -fp-model source:
  OPT_test   := $(OPT_normal)
  # provide useful info at exceptions (use for debug):
  OPT_debug  := -O0 -traceback -check bounds -check uninit -ftrapuv
  # provides many warnings:
  OPT_warn   := -warn

  mode = normal
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

  MPIIFORT := $(shell command -v mpiifort 2> /dev/null)
  ifdef MPIFC
    COMP_PAR = $(MPIFC)
  endif
  ifndef MPIFC
    ifdef MPIIFORT
      COMP_PAR = mpiifort
    endif
    ifndef MPIIFORT
      COMP_PAR = mpifort
    endif
  endif

  LINK_PAR = $(COMP_PAR)
  OPT_PAR  = $(OPT_SER)

  LLIB =

endif

#------------------
#    LOCAL_GFORTRAN
#------------------

ifeq ($(ARCH), LOCAL_GFORTRAN)

  COMP_SER = gfortran
  LINK_SER = $(COMP_SER)
  OPT_SER  = -O3 -ffree-form -ffree-line-length-none

  OPT_format := -ffree-form -ffree-line-length-none
  OPT_normal := $(OPT_format) -O2
  # quick compilation:
  OPT_quick  := $(OPT_format) -O0
  OPT_test   := $(OPT_normal)
  # provide useful info at exceptions (use for debug):
  OPT_debug  := $(OPT_format) -fimplicit-none -fbacktrace -g -fcheck=all -ffpe-trap=zero,overflow -finit-real=nan
  # provides many warnings:
  OPT_warn   := $(OPT_format) -O0 -Wall -Wno-character-truncation -Wno-conversion -Wno-unused-function
  # time the execution
  OPT_gprof  := $(OPT_format) -O0 -pg

  mode = normal
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
  ifeq ($(mode),gprof)
    OPT_SER = $(OPT_gprof)
  endif
  LLIB =

  GCCVERSIONLTEQ4 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \<= 4)
  ifeq "$(GCCVERSIONLTEQ4)" "1"
    OPT_SER += -D_NOIEEE_
  endif
  ifdef MPIFC
    COMP_PAR = $(MPIFC)
  endif
  ifndef MPIFC
    COMP_PAR = mpifort
  endif
  LINK_PAR = $(COMP_PAR)
  OPT_PAR  = $(OPT_SER)


endif

#-----------------
#     ALARIK_INTEL
#-----------------

ifeq ($(ARCH), ALARIK_INTEL)

  COMP_SER = ifort
  LINK_SER = $(COMP_SER)
  OPT_SER  = -O3

  COMP_PAR = mpiifort
  LINK_PAR = $(COMP_PAR)
  OPT_PAR  = $(OPT_SER)

  LLIB =

# for Intel compiler on Alarik we need static linking for FFTW3 to work
#    LINK_SER += -Bstatic
#    LINK_PAR += -Bstatic

endif

#--------------------
#     ALARIK_GFORTRAN  (Keep ARCH = ALARIK_INTEL for NOW)
#--------------------

ifeq ($(ARCH), ALARIK_GFORTRAN)

  COMP_SER = gfortran
  LINK_SER = $(COMP_SER)
  OPT_SER  = -O3 -ffree-form -ffree-line-length-none  -g

  COMP_PAR = mpifort
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
  ifeq ($(mode),normal)
    DEFINES_BASE += -D_NORMAL_
  endif
  ifeq ($(mode),debug)
    DEFINES_BASE += -D_DEBUG_
  endif
  ifeq ($(mode),warn)
    DEFINES_BASE += -D_WARN_
  endif
  ifeq ($(mode),quick)
    DEFINES_BASE += -D_QUICK_
  endif
  ifeq ($(mode),gprof)
    DEFINES_BASE += -D_GPROF_
  endif
endif


DEFINES_SER  = $(DEFINES_BASE)
DEFINES_PAR  = $(DEFINES_BASE) -D_PAR_

#-----------------------------------------------------------------
#                          Compiler flags
#-----------------------------------------------------------------

FLAG_SER  = $(DEFINES_SER) $(OPT_SER)
FLAG_PAR  = $(DEFINES_PAR) $(OPT_PAR)

FLAG_SER += -I$(FFTW_PATH)
FLAG_PAR += -I$(FFTW_PATH)

mollib = $(abspath molsim.lib)

#-----------------------------------------------------------------
#                          Common libraries
#-----------------------------------------------------------------

ifeq ($(FFTWLIB), none)
  FFTWLIB=$(FFTW_PATH)/lib
endif

ifeq ($(ARCH), LOCAL_INTEL)
  LLIB += -L$(FFTWLIB) -lfftw3 -lm
endif
ifeq ($(ARCH), LOCAL_GFORTRAN)
  LLIB += -L$(FFTWLIB) -lfftw3 -lm
endif
ifeq ($(ARCH), ALARIK_INTEL)
  LLIB += $(FFTW3_HOME)/lib/libfftw3.so  #linker issues with ifort on Alarik => workaround
endif


#-----------------------------------------------------------------

OBJS1=mol.o            \
      particle.o       \
      potential.o      \
      coordinate.o     \
      md.o             \
      mc.o             \
      bd.o             \
      nlist.o          \
      energy.o         \
      denergy.o        \
      dump.o           \
      group.o          \
      static.o         \
      dynamic.o        \
      image.o          \
      statistics.o     \
      mixed.o          \
      molaux.o         \
      mollib.o         \
      moluser.o        \
      sso.o            \
      mesh.o           \
      celllist.o

OBJS2=$(OBJS1:.o=_par.o)
OBJS2+=parallel_par.o

OBJS3=mol.o            \
      statistics.o     \
      molaux.o         \
      mollib.o         \
      mesh.o

.PHONY: default def all molsim_ser molsim_par ser par moldyn ver install uninstall clean tidy testmesh

default: def
def:
	@echo "Specify any one: of all, ser, par, molsim_ser, molsim_par, moldyn"
all: ser par
molsim_ser: $(BINDIR)/molsim_ser.exe
molsim_par: $(BINDIR)/molsim_par.exe
ser: $(BINDIR)/molsim_ser.exe
par: $(BINDIR)/molsim_par.exe
moldyn: $(BINDIR)/moldyn.exe
ver:
	@[[ -z  "$(VER)"  ]] && echo No version name has been specified. || echo Version is specified as \"$(VER)\"
	@echo bindir is \"$(BINDIR)\"

install:
	cp molsim_ser $(HOME)/bin/molsim_ser$(VER); sed -i "s|^bin=.*|bin=$(BINDIR)|g" $(HOME)/bin/molsim_ser$(VER)
	cp molsim_par $(HOME)/bin/molsim_par$(VER); sed -i "s|^bin=.*|bin=$(BINDIR)|g" $(HOME)/bin/molsim_par$(VER)

uninstall:
	rm -f $(HOME)/bin/molsim_par$(VER) $(HOME)/bin/molsim_ser$(VER)

$(BINDIR)/molsim_ser.exe: molsim.o $(OBJS1)
	mkdir -p $(BINDIR)
	$(LINK_SER) $(FLAG_SER) -o $@ $^ $(LLIB)

$(BINDIR)/molsim_par.exe: molsim_par.o $(OBJS2)
	mkdir -p $(BINDIR)
	$(LINK_PAR) $(FLAG_PAR) -o $@ $^ $(LLIB)

$(BINDIR)/moldyn.exe:     moldyn.o $(OBJS3)
	$(COMP_SER) $(FLAG_SER) -o moldyn.exe moldyn.o $(OBJS3) $(LLIB)
	mkdir -p $(BINDIR); mv moldyn.exe $(BINDIR)/moldyn.exe; cp moldyn $(HOME)/bin/moldyn$(VER)

testmesh:   testmesh.o mesh.o
	$(COMP_SER) $(FLAG_SER) -o testmesh.exe testmesh.o mesh.o
	./testmesh.exe

%.o: %.F90
	$(COMP_SER) $(FLAG_SER) -o $@ -c $<
%_par.o: %.F90
	$(COMP_PAR) $(FLAG_PAR) -o $@ -c $<

mol.o: FLAG_SER += -DFLIBMACRO=\"$(mollib)\"
mol_par.o: FLAG_PAR += -DFLIBMACRO=\"$(mollib)\"

bd.o:             mol.o md.o
bd_par.o:         mol_par.o md_par.o
celllist.o:       mol.o
celllist_par.o:   mol_par.o
coordinate.o:     mol.o mollib.o
coordinate_par.o: mol_par.o mollib_par.o
denergy.o:        mol.o energy.o mollib.o celllist.o
denergy_par.o:    mol_par.o energy_par.o mollib_par.o celllist_par.o
dump.o:           mol.o
dump_par.o:       mol_par.o
dynamic.o:        mol.o statistics.o mollib.o
dynamic_par.o:    mol_par.o statistics_par.o mollib_par.o
energy.o:         mol.o mollib.o celllist.o
energy_par.o:     mol_par.o mollib_par.o celllist_par.o
group.o:          mol.o
group_par.o:      mol_par.o
image.o:          mol.o mollib.o
image_par.o:      mol_par.o mollib_par.o
mc.o:             mol.o nlist.o mollib.o celllist.o
mc_par.o:         mol_par.o nlist_par.o mollib_par.o celllist_par.o
md.o:             mol.o
md_par.o:         mol_par.o
mesh.o:
mesh_par.o:
mixed.o:          mol.o
mixed_par.o:      mol_par.o
mol.o:            statistics.o mesh.o mollib.o
mol_par.o:        statistics_par.o mesh_par.o mollib_par.o
molaux.o:         mol.o mollib.o
molaux_par.o:     mol_par.o mollib_par.o
moldyn.o:         statistics.o
moldyn_par.o:     statistics_par.o
mollib.o:
mollib_par.o:
molsim.o:         mol.o statistics.o mollib.o
molsim_par.o:     mol_par.o statistics_par.o mollib_par.o
moluser.o:        mol.o statistics.o potential.o mollib.o
moluser_par.o:    mol_par.o potential_par.o statistics_par.o mollib_par.o
nlist.o:          mol.o celllist.o
nlist_par.o:      mol_par.o celllist_par.o
parallel_par.o:   mol_par.o
particle.o:       mol.o mollib.o
particle_par.o:   mol_par.o mollib_par.o
potential.o:      mol.o
potential_par.o:  mol_par.o
sso.o:            mc.o mol.o
sso_par.o:        mc_par.o mol_par.o
static.o:         mol.o statistics.o mollib.o potential.o
static_par.o:     mol_par.o statistics_par.o mollib_par.o potential_par.o
statistics.o:     mollib.o
statistics_par.o: mollib_par.o
testmesh.o:       mesh.o
testmesh_par.o:   mesh_par.o

#*****************************************************************
#*                         Tidy                                  *
#*****************************************************************

tidy:
	rm -f *.o *.mod

#*****************************************************************
#*                         Clean                                 *
#*****************************************************************

clean:
	rm -f *.o *.mod *.i90 $(BINDIR)/*.exe *__genmod*
