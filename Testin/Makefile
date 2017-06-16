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

todifffile := todiff.txt #TODO: change to real todiff file
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
.PHONY: save declarestable clean cleansave cleanout

########
# diff #
########

all: diff.out

diff.out:	$(difffiles)
	cat $^ > $@

########
# diff #
########

$(OutDir)/%.diff: $(OutDir)/% $(SaveDir)/%
	$(diffscript) $^ > $@

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

cleanall: clean cleansave

################
# special runs #
################

%.ana.out: %.pos %.ana.in
	for i in pos ori liv anv for tor group utot; do rm -f $*.ana.$$i; ln -s $$(pwd)/$*.$$i $*.ana.$$i ; done
	cd $(@D); $(molsim) $(basename $(@F))
	for i in pos ori liv anv for tor group utot; do rm -f $*.ana.$$i; done

%/dip.md.ew.out: %/dip.md.ew.cnf %/dip.md.ew.in
	cd $(@D); $(molsim) $(basename $(@F))

%/dip.md.ew.cnf:
	cp $(InDir)/dip.md.ew.cnf.save $@

%/daniel.out: %/daniel.cnf %/daniel.in
	cd $(@D); $(molsim) $(basename $(@F))

%/daniel.cnf:
	cp $(InDir)/daniel.cnf.save $@

%/ads_20-mer:	%/p.bd.ads2.out
	@echo "ads_20-mer is created"

%/p.bd.ads2.out %/p.bd.ads2.user: %/p.bd.ads2.cnf %/p.bd.ads2.in
	cd $(@D); $(molsim) $(basename $(@F))
	rm $<

%/p.bd.ads2.cnf: %/p.bd.ads1.cnf
	cp $< $@ ; cp $*/p.bd.ads1.user $*/p.bd.ads2.user

$(OutDir)/pb.bd.asd2.user.diff: $(OutDir)/p.bd.ads2.user $(OutDir)/p.bd.ads.user
	$(diffscript) $^ > $@

%/sb_xz.out: %/sb_xz.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "sb_xz exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

################
# molsim files #
################

%.out %.wrl %.pos %.cnf %.user: %.in
	cd $(@D); $(molsim) $(*F)

$(OutDir)/%.in: $(InDir)/%.in
	$(mkdirp) $(@D)
	cp $< $@

$(SaveDir)/%.in: $(InDir)/%.in
	$(mkdirp) $(@D)
	cp $< $@
