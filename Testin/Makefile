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
.PHONY: save declarestable clean cleansave cleanout checkcurrent molsim_exe diff out

########
# diff #
########

all: diff.out

diff: diff.out

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

$(savefiles): | checkcurrent

checkcurrent: $(stable) $(current)
	@cmp --silent $(stable) $(current) || { echo "the current version is not the stable one"; exit 1 ; }

$(current): $(sourcefiles)
	@echo $(sourcefiles) > $@
	@cat $(sourcefiles) | md5sum >> $@

declarestable:
	@echo $(sourcefiles)  > $(stable)
	@cat $(sourcefiles) | md5sum >> $(stable)

#############
# out files #
#############

out: $(outfiles)

##########
# molsim #
##########

$(outfiles) $(savefiles): | molsim_exe

molsim_exe:
	@echo recompiling molsim_ser
	cd ../Src/ && make clean && make ser mode=test && cd ../Testin


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
	@echo ""

%/p.bd.ads2.out %/p.bd.ads2.user: %/p.bd.ads2.cnf %/p.bd.ads2.in
	cd $(@D); $(molsim) $(basename $(@F))
	rm $<

%/p.bd.ads2.cnf: %/p.bd.ads1.cnf
	cp $< $@ ; cp $*/p.bd.ads1.user $*/p.bd.ads2.user

$(OutDir)/pb.bd.asd2.user.diff: $(OutDir)/p.bd.ads2.user $(OutDir)/p.bd.ads.user
	$(diffscript) $^ > $@

%/sb_xz.out: %/sb_xz.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

%/sb_x.out: %/sb_x.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

%/sb_xyz.out: %/sb_xyz.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

%/sb_y.out: %/sb_y.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

%/sb_z.out: %/sb_z.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

%/sb_z_y090.out: %/sb_z_y090.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

%/sb_z_y180.out: %/sb_z_y180.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

%/sb_z_y030.out: %/sb_z_y030.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

%/sb_z_y045.out: %/sb_z_y045.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

%/sb_z_xyz.out: %/sb_z_xyz.in
	cd $(@D); $(molsim) $(basename $(@F)); if [ $$? = 1 ] ; then echo "exited with exit code 1, as expected"; else echo "sb_xz did not finish as exptected"; exit 1; fi

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
