CONFFILES=Src/make.arch
all: ser par install

.PHONY: ser par doc install clean

doc:
	+(cd Doc/; $(MAKE) doc)

install:
	+(cd Src/; $(MAKE) install)

ser: $(CONFFILES)
	+(cd Src/; $(MAKE) ser)

par: $(CONFFILES)
	+(cd Src/; $(MAKE) par)

$(CONFFILES): configure.sh
	@echo "you need to run the configure script first! (./configure.sh)" && exit 1

clean:
	+(cd Src/; $(MAKE) clean)
	+(cd Doc/; $(MAKE) clean)
