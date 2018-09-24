.SUFFIXES:
.INTERMEDIATE: latex/refman.pdf latex/refman.tex
.PHONY: documentation doc clean
documentation: doc
doc: ../documentation.pdf

../documentation.pdf: latex/refman.pdf
	cp $< $@

latex/refman.pdf: latex/refman.tex
	cd latex && make pdf

latex/refman.tex: doxy header.tex footer.tex *.md ../Src/*.F90
	doxygen $<

clean:
	rm -rf latex
