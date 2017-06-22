#*****************************************************************
#*       Makefile for the MOLSIM Tests                           *
#*****************************************************************
# This Makefile should be used to generate output form the current stable version (make save) and then to compare the results of a modified version of MOLSIM to the output of the stable version (make diff)

###########
# Scripts #
###########

# the script to diff two output files, will crop out trivial differences like version numbers and time stamps
diffscript = scripts/diff1.sh
# this script will run MOLSIM. When creating something in the save directory, it is checked, that the current version is the stable one
molsim = scripts/molsim_save.sh

#mechanism to get the correctly named MOLSIM script
versionfile = ../version.conf
molsim_bin = ../Bin/molsim_ser.exe

#commands to make directories and delete stuff
rm = rm -rf
mkdirp = mkdir -p

########
# Dirs #
########

#when changing these variables, modify them also in scripts/molsim_save.sh
SaveDir = save
OutDir = out
InDir = in

#########
# files #
#########

#which files are to be diffed is stored in the following files:
todifffile := todiff.txt

#get the input, output and save filenames
todiff := $(shell cat $(todifffile))
infiles := $(wildcard $(InDir)/*.in)

outfiles = $(addprefix $(OutDir)/,$(todiff))
outinfiles = $(subst $(InDir), $(OutDir) ,$(infiles))

savefiles = $(addprefix $(SaveDir)/,$(todiff))
saveinfiles = $(subst $(InDir), $(SaveDir) ,$(infiles))

#store the difference between the save and out file in %.diff
difffiles = $(addsuffix .diff,$(outfiles))

#all source files which have an effect on the final program
sourcefiles_un := $(wildcard ../Src/*.F90)
sourcefiles_un += ../Src/makefile ../configure.sh ../Src/molsim_ser
sourcefiles = $(sort $(sourcefiles_un))

intermediateSpecialFiles = $(SaveDir)/p.bd.ads1.cnf $(OutDir)/p.bd.ads1.cnf

#files to define the stable and the current version
# also adapt these variables in script/molsim_save.sh
stable=stable.md5sum
current=.current.md5sum

#file where all diffs are collected
diffout = diff.out

# these files are not necessarily needed, when all prerequisites are met, but we do not want to delete them automatically:
.SECONDARY: $(outfiles) $(outinfiles) $(difffiles)
#we always need to make sure molsim is compiled, and store when we compiled it in the timestamp of $(compiled):
.PRECIOUS: $(current) $(savefiles) $(intermediateSpecialFiles)
#when molsim exists with an error, do not keep the faulty files
.DELETE_ON_ERROR:
#these are only used to interact with the makefile
.PHONY: all diff save out declarestable clean cleansave cleanout
# use diff to get the differences between the stable and the current version
# use save to compile the save directory. Do this when the current stable version is represented in the source files (the master branch)
# use out to generate the output files with the current version
# use decarestable to declare the current version as a stable one
# use clean to delete the out directory, the $(current) and the $(diffout) file
# use cleansave to delete the save directory  the $(current) and the $(diffout) file
# use cleanall to run clean and cleansave

########
# diff #
########

all: $(diffout)

diff: $(diffout)

$(diffout):	$(difffiles)
	cat $^ > $@

$(OutDir)/%.diff: $(OutDir)/% $(SaveDir)/%
	$(diffscript) $^ > $@

##############
# save files #
##############

save: $(savefiles)

#when the $(stable) file is not existent, then we do not know which version is stable, and we cannot create the save dir
$(savefiles): $(stable)

declarestable: $(current)
	cp $(current) $(stable)

#############
# out files #
#############

out: $(outfiles)

$(outfiles): $(sourcefiles)

###################
# current version #
###################

# all files need to have the current version defined
$(outfiles) $(savefiles): $(current)

# when the source files are newer than the molsim_ser script (or the script is absent), than mark the script as new, which triggers the $(compiled) target.
$(molsim_bin): $(sourcefiles)
	@touch $(molsim_bin)

# when $(current) does not exist, it is assumed that a recompilation is needed. From then on, molsim_bin is older than $(current) except for the case where the source files are newer than the script (which triggers a touch of the script) or MOLSIM was recompiled. Under these circumstances, the current state of the compilation is unknown, and therefore we need to recompile.
#how to make the current file
$(current): $(molsim_bin) $(sourcefiles)
	+cd ../Src/ && make clean && make ser mode=test && cd ../Testin
	echo $(sourcefiles)  > $(current)
	cat $(sourcefiles) | md5sum >> $(current)
	$(molsim_bin) -v >> $(current)


#########
# clean #
#########

#various clean commands
#by default only clean the out dir
clean:
	$(rm) $(OutDir) $(current) $(diffout)

cleansave:
	$(rm) $(SaveDir) $(current) $(diffout)

cleanall: clean cleansave

################
# special runs #
################

#define some rules on how to prepare some special files
%.ana.out: %.pos %.ana.in
	for i in pos ori liv anv for tor group utot; do rm -f $*.ana.$$i; ln -s $$(pwd)/$*.$$i $*.ana.$$i ; done
	$(molsim) $@
	for i in pos ori liv anv for tor group utot; do rm -f $*.ana.$$i; done

%/dip.md.ew.out: %/dip.md.ew.cnf %/dip.md.ew.in
	$(molsim) $@

%/dip.md.ew.cnf:
	cp $(InDir)/dip.md.ew.cnf.save $@

%/daniel.out: %/daniel.cnf %/daniel.in
	$(molsim) $@

%/daniel.cnf:
	cp $(InDir)/daniel.cnf.save $@

$(OutDir)/p.bd.ads2.user.diff: $(OutDir)/p.bd.ads2.user $(OutDir)/p.bd.ads.user
	$(diffscript) $^ > $@

%/ads_20-mer: %/p.bd.ads.out
	touch $@

%/p.bd.ads2.out: %/p.bd.ads2.in %/p.bd.ads1.in
	$(molsim) $*/p.bd.ads1.out
	mv $*/p.bd.ads1.cnf $*/p.bd.ads2.cnf
	mv $*/p.bd.ads1.user $*/p.bd.ads2.user
	$(molsim) $*/p.bd.ads2.out

################
# MOLSIM files #
################

# default rule on how to run MOLSIM
%.out %.wrl %.pos %.cnf %.user: %.in $(current)
	$(molsim) $@

$(OutDir)/%.in: $(InDir)/%.in
	$(mkdirp) $(@D)
	cp $< $@

$(SaveDir)/%.in: $(InDir)/%.in
	$(mkdirp) $(@D)
	cp $< $@
