#*****************************************************************
#*       Makefile for the MOLSIM Tests                           *
#*****************************************************************

###########
# Scripts #
###########

diffscript = scripts/diff1.sh

versionfile = ../version.conf
ifneq ("$(wildcard $(versionfile))","")
version := `cat $(versionfile)`
else
version :=
endif
ifeq ($(version),)
molsim = "$$HOME/bin/molsim_ser.$(version)"
else
molsim = "$$HOME/bin/molsim_ser"
endif

rm = rm -rf
mkdirp = mkdir -p

########
# Dirs #
########

SaveDir = Save
OutDir = Out
InDir = In

#########
# files #
#########

todifffile := todiff_test.txt #TODO: change to real todiff file
todiff := $(shell cat $(todifffile))
infiles := $(wildcard $(InDir)/*.in)

outfiles = $(addprefix $(OutDir)/,$(todiff))
outinfiles = $(subst $(InDir), $(OutDir) ,$(infiles))

savefiles = $(addprefix $(SaveDir)/,$(todiff))
saveinfiles = $(subst $(InDir), $(SaveDir) ,$(infiles))

difffiles = $(addsuffix .diff,$(outfiles))

sourcefiles := $(wildcard ../Src/*.F90)
sourcefiles += ../Src/makefile ../configure.sh

stable = stable.md5sum
current = .current.md5sum

.PRECIOUS: $(outfiles) $(outinfiles) $(savefiles) $(difffiles)
.PHONY: save declarestable

########
# diff #
########

all: diff.out

diff.out:	$(difffiles)
	cat $^ > $@

########
# diff #
########

$(OutDir)/%.diff: $(OutDir)/%
	$(diffscript) $< $(SaveDir)/$* > $@

################
# molsim files #
################

%.out %.wrl: %.in
	cd $(@D); $(molsim) $(*F)

$(OutDir)/%.in $(SaveDir)/%.in: $(InDir)/%.in
	$(mkdirp) $(@D)
	cp $< $@

##############
# save files #
##############

save: $(savefiles)
	$(rm) .StableVer

$(savefiles): .StableVer

.StableVer: $(stable) $(current)
	@cmp --silent $(stable) $(current) || { echo "the current version is not the stable one"; exit 1 ; }
	echo "currently_stable" >  $@

$(current): $(sourcefiles)
	@echo $(sourcefiles) > $@
	@cat $(sourcefiles) | md5sum >> $@

declarestable:
	@echo $(sourcefiles)  > $(stable)
	@cat $(sourcefiles) | md5sum >> $(stable)

#########
# clean #
#########

clean:
	$(rm) $(OutDir)

cleansave:
	$(rm) $(SaveDir)
