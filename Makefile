#
#  Makefile by Scott Prahl, Mar 2024
#

export VERSION = 3-16-3

#Base directory for installation
DESTDIR=/usr/local

#typical install hierarchy
BIN_INSTALL=$(DESTDIR)/bin
LIB_INSTALL=$(DESTDIR)/lib
INC_INSTALL=$(DESTDIR)/include

#Library extension
LIB_EXT = .dylib      # macOS shared lib
#LIB_EXT = .so        # linux shared lib
#LIB_EXT = .a         # static lib

# default executable name
IAD_EXECUTABLE = ./iad
AD_EXECUTABLE = ./ad

CFLAGS = -dynamic -fno-common -Wall -ansi  #First two flags needed on macOS to build .dylib

export IAD_OBJ = src/iad_util.o  src/iad_calc.o src/iad_find.o src/iad_pub.o  src/iad_io.o

export AD_OBJ  = src/nr_zbrak.o  src/ad_bound.o src/ad_doubl.o src/ad_frsnl.o src/ad_globl.o \
                 src/ad_matrx.o  src/ad_phase.o src/ad_prime.o src/ad_radau.o src/ad_start.o \
                 src/ad_cone.o   src/ad_layers.o

export NR_OBJ  = src/nr_amoeb.o  src/nr_amotr.o src/nr_brent.o src/nr_gaulg.o src/nr_mnbrk.o \
                src/nr_rtsaf.o  src/nr_util.o  src/nr_hj.o

MAIN = Makefile INSTALL.md README.md License

DOCS =  CHANGELOG.rst         docs/ToDo.md               docs/manual.tex      \
        docs/ad_src.pdf       docs/iad_src.pdf           docs/manual.pdf      \
        docs/ch3RTcorr.pdf    docs/ch3spheremeas.pdf     docs/ch3spheresR.pdf \
        docs/ch3spheresT.pdf  docs/ch3Doublespheres.pdf  docs/colltrans.pdf   \
        docs/lightloss.pdf    docs/niek_graph.pdf        docs/glass_slide.pdf \
        docs/cmdexe.png       docs/valid.png             docs/dual8.png       \
        docs/dual90.png       docs/iad.bib

TEST =  test/Makefile       test/basic_A.rxt    test/basic_B.rxt      test/basic_C.rxt    test/basic_D.rxt    \
        test/double.rxt     test/example2.rxt   test/il_A.rxt         test/il_B.rxt       test/il_C.rxt       \
        test/ink_A.rxt      test/ink_B.rxt      test/ink_C.rxt        test/kenlee_A.rxt   test/kenlee_B.rxt   \
        test/kenlee_C.rxt   test/newton.rxt     test/royston2.rxt     test/royston3_A.rxt test/royston3_B.rxt \
        test/royston3_C.rxt test/royston3_D.rxt test/royston3_E.rxt   test/royston9_A.rxt test/royston9_B.rxt \
        test/royston9_C.rxt test/royston9_D.rxt test/royston1.rxt     test/sample_A.rxt   test/sample_B.rxt   \
        test/sample_C.rxt   test/sample_D.rxt   test/sample_E.rxt     test/sample_F.rxt   test/sample_G.rxt   \
        test/sevick_A.rxt   test/sevick_B.rxt   test/terse_A.rxt      test/terse_B.rxt    test/tio2_vis.rxt   \
        test/uterus.rxt     test/valid.bat      test/vio_A.rxt        test/vio_B.rxt      test/x_bad_data.rxt \
        test/ville1.rxt     test/fairway_A.rxt  test/fairway_B.rxt    test/fairway_C.rxt  test/fairway_D.rxt  \
        test/fairway_E.rxt  test/basic_E.rxt    test/combo_0.rxt

export WSRC =  src/ad.w     src/ad_frsnl.w       src/ad_prime.w        src/iad_io.w          \
        src/ad_globl.w      src/ad_radau.w       src/iad_main.w                              \
        src/ad_bound.w      src/ad_layers.w      src/ad_start.w                              \
        src/ad_chapter.w    src/ad_layers_test.w src/iad.w             src/iad_pub.w         \
        src/ad_cone.w       src/ad_main.w        src/iad_calc.w        src/iad_type.w        \
        src/ad_cone_test.w  src/ad_matrx.w       src/iad_util.w        src/ad_oblique_test.w \
        src/ad_doubl.w      src/ad_phase.w       src/iad_find.w        src/iad_agrid.w       \
        src/mc_lost.w

export NRSRC = src/nr_amoeb.c src/nr_amotr.h     src/nr_gaulg.c        src/nr_mnbrk.h  \
        src/nr_util.c       src/nr_util.h        src/nr_zbrak.c        src/nr_zbrak.h  \
        src/nr_amoeb.h      src/nr_brent.c       src/nr_gaulg.h        src/nr_rtsaf.c  \
        src/nr_amotr.c      src/nr_brent.h       src/nr_mnbrk.c        src/nr_rtsaf.h  \
        src/nr_hj.c         src/nr_hj.h          src/nr_zbrent.h       src/nr_zbrent.c \
        src/mc_lost_main.c  src/version.h        src/mc_test.c

export CSRC  = src/ad_frsnl.c      src/ad_globl.c       src/ad_matrx.c        src/ad_start.c        \
        src/iad_main.c      src/ad_doubl.c       src/iad_util.c        src/ad_radau.c        \
        src/ad_prime.c      src/iad_find.c       src/ad_phase.c        src/ad_bound.c        \
        src/ad_layers.c     src/version.c        src/iad_io.c          src/ad_chapter.c      \
        src/iad_calc.c      src/iad_pub.c        src/ad_cone.c         src/ad_oblique_test.c \
        src/ad_cone_test.c  src/ad_layers_test.c src/iad_agrid.c       src/mc_lost.c

export HSRC  = src/ad_bound.h      src/ad_globl.h       src/ad_phase.h        src/ad_start.h   src/iad_io.h   \
        src/ad_doubl.h                           src/ad_prime.h        src/iad_calc.h   src/iad_util.h \
        src/ad_frsnl.h      src/ad_matrx.h       src/ad_radau.h        src/iad_find.h   src/iad_pub.h  \
        src/ad_cone_ez.h    src/ad_layers.h      src/ad_cone.h         src/iad_type.h   src/iad_agrid.h \
        src/mc_lost.h

OSRC  = src/system.bux src/ad.bux src/iad.bux src/cobweb.pl src/version.pl src/Makefile src/toDOS.pl

all :  ad iad

lib :
	cd src ; make libiad.h
	cp src/libiad.h .
	cd src ; make libiad$(LIB_EXT)
	cp src/libiad$(LIB_EXT) .

install: ad iad
	mkdir -p $(BIN_INSTALL)
	rm -f $(BIN_INSTALL)/ad
	cp ad  $(BIN_INSTALL)
	rm -f $(BIN_INSTALL)/iad
	cp -f iad $(BIN_INSTALL)

install-lib: lib libiad$(LIB_EXT) libiad.h
	mkdir -p $(LIB_INSTALL)
	mkdir -p $(INC_INSTALL)
	cp libiad.h $(INC_INSTALL)
	cp libiad$(LIB_EXT) $(LIB_INSTALL)

dists: dist windist

dist:
	touch src/version.h
	cd src && ./version.pl
	make docs
	make clean
	make
	make tidy
	mkdir -p    iad-$(VERSION)
	mkdir -p    iad-$(VERSION)/docs
	mkdir -p    iad-$(VERSION)/test
	mkdir -p    iad-$(VERSION)/src
	ln $(MAIN)  iad-$(VERSION)
	ln $(DOCS)  iad-$(VERSION)/docs
	ln $(TEST)  iad-$(VERSION)/test
	ln $(WSRC)  iad-$(VERSION)/src
	ln $(HSRC)  iad-$(VERSION)/src
	ln $(CSRC)  iad-$(VERSION)/src
	ln $(NRSRC) iad-$(VERSION)/src
	ln $(OSRC)  iad-$(VERSION)/src
	zip -r archives/iad-$(VERSION) iad-$(VERSION)
	rm -rf iad-$(VERSION)
	cp archives/iad-$(VERSION).zip archives/iad-latest.zip

windist: ad.exe iad.exe libiad.dll
	make docs
	make tidy
	mkdir -p      iad-win-$(VERSION)
	mkdir -p      iad-win-$(VERSION)/docs
	mkdir -p      iad-win-$(VERSION)/test
	mkdir -p      iad-win-$(VERSION)/src
	ln ad.exe     iad-win-$(VERSION)
	ln iad.exe    iad-win-$(VERSION)
	ln libiad.dll iad-win-$(VERSION)
	ln $(MAIN)    iad-win-$(VERSION)
	ln $(DOCS)    iad-win-$(VERSION)/docs
	ln $(TEST)    iad-win-$(VERSION)/test
	ln $(WSRC)    iad-win-$(VERSION)/src
	ln $(HSRC)    iad-win-$(VERSION)/src
	ln $(CSRC)    iad-win-$(VERSION)/src
	ln $(NRSRC)   iad-win-$(VERSION)/src
	ln $(OSRC)    iad-win-$(VERSION)/src
	src/toDOS.pl  iad-win-$(VERSION)/src/*.c
	src/toDOS.pl  iad-win-$(VERSION)/src/*.h
	src/toDOS.pl  iad-win-$(VERSION)/test/*.rxt
	src/toDOS.pl  iad-win-$(VERSION)/test/valid.bat
	rm iad-win-$(VERSION)/src/*.bak
	rm iad-win-$(VERSION)/test/*.bak
	zip -r archives/iad-win-$(VERSION) iad-win-$(VERSION)
	rm -rf iad-win-$(VERSION)
	cp archives/iad-win-$(VERSION).zip archives/iad-win-latest.zip

ad: $(WSRC) $(NRSRC)
	cd src ; make ad
	cp src/ad ad

iad: $(WSRC) $(NRSRC)
	cd src ; make iad
	cp src/iad iad

ad.exe: $(WSRC) $(NRSRC)
	cd src ; make clean
	cd src ; make CC="x86_64-w64-mingw32-gcc" ad
	mv src/ad.exe ad.exe
	cd src ; make clean

iad.exe: $(WSRC) $(NRSRC)
	cd src ; make clean
	cd src ; make CC="x86_64-w64-mingw32-gcc" iad
	mv src/iad.exe iad.exe
	cd src ; make clean

libiad.dll: $(WSRC) $(NRSRC)
	cd src ; make clean
	cd src ; make CC="x86_64-w64-mingw32-gcc" libiad.dll
	mv src/libiad.dll .
	cd src ; make clean

docs:
	perl -pi -e 's/\\def\\adversion.*/\\def\\adversion{$(VERSION)}/' src/ad.w
	cd src ; make ad_doc
	perl -pi -e 's/\\def\\iadversion.*/\\def\\iadversion{$(VERSION)}/' src/iad.w
	cd src ; make iad_doc
	cd docs ; pdflatex manual
	cd docs ; bibtex manual
	cd docs ; pdflatex manual
	cd docs ; pdflatex manual

clean:
	rm -f ad iad *.pdf libiad.a libiad.so libiad.dylib libiad.h src/lib_iad.h src/lib_ad.h
	rm -f src/*.o tests/*.abg
	rm -f src/*.aux src/*.dvi src/*.idx src/*.ref src/*.sref src/*.tex src/*.toc src/*.log src/*.scn
	rm -f iad.exe ad.exe src/iad.exe src/ad.exe
	rm -f libiad.dll src/libiad.dll
	rm -f src/oblique_test src/mc_test src/cone_test src/layer_test src/mc_lost
	rm -rf tests/.jupyter tests/.ipynb_checkpoints .jupyter .ipynb_checkpoints

clean-generated-dry-run:
	tools/clean_generated.py

clean-generated:
	tools/clean_generated.py --execute

clean-generated-all-dry-run:
	tools/clean_generated.py --include-tracked-generated

clean-generated-all:
	tools/clean_generated.py --execute --include-tracked-generated

clean-generated-archives-dry-run:
	tools/clean_generated.py --include-archives

clean-generated-archives:
	tools/clean_generated.py --execute --include-archives

scratch-build: clean-generated-all
	$(MAKE) tidy
	$(MAKE) docs
	$(MAKE) all

realclean:
	make clean
	cd src ; make realclean
	rm -f ad iad libiad.h libiad$(LIB_EXT)
	rm -f docs/manual.pdf docs/manual.bbl docs/manual.blg docs/manual.out docs/manual.toc
	rm -f docs/manual.aux docs/manual.log
	cd tests ; make clean
	rm -f $(CSRC) $(HSRC)

tidy:
	$(MAKE) -C src tidy-src

mc_lost:
	cd src; make mc_lost
	cp src/mc_lost .

mc_test:
	cd src ; make mc_test
	src/mc_test

mc_lost_test: mc_lost
	cd src ; make mc_lost
	src/mc_lost
	src/mc_lost -P 0
	src/mc_lost -i 0 -n 1.0
	src/mc_lost -i 0 -n 2
	src/mc_lost -n 1.4
	src/mc_lost -n 1.4 -N 1.5
	src/mc_lost -n 1.4 -N 1.5 -g 0.9

layer_test:
	cd src ; make layer_test
	src/mc_test

cone_test:
	cd src ; make cone_test
	src/cone_test

oblique_test:
	cd src ; make oblique_test
	src/oblique_test

executables:
	cd src; make ad
	cd src; make iad
	cd src; make mc_test
	cd src; make layer_test
	cd src; make oblique_test
	cd src; make cone_test
	cd src; make mc_lost

veryshorttest: iad
	IAD_EXECUTABLE=$(IAD_EXECUTABLE) tests/cli/run_cli_tests.sh veryshort

test shorttest: iad ad
	IAD_EXECUTABLE=$(IAD_EXECUTABLE) AD_EXECUTABLE=$(AD_EXECUTABLE) tests/cli/run_cli_tests.sh basic

longtest: iad
	IAD_EXECUTABLE=$(IAD_EXECUTABLE) tests/cli/run_cli_tests.sh batch tests/rxt nomc

test0: iad
	IAD_EXECUTABLE=$(IAD_EXECUTABLE) tests/cli/run_cli_tests.sh batch tests/rxt/0_sphere

test1: iad
	IAD_EXECUTABLE=$(IAD_EXECUTABLE) tests/cli/run_cli_tests.sh batch tests/rxt/1_sphere

test2: iad
	IAD_EXECUTABLE=$(IAD_EXECUTABLE) tests/cli/run_cli_tests.sh batch tests/rxt/2_sphere

test1_nomc: iad
	IAD_EXECUTABLE=$(IAD_EXECUTABLE) tests/cli/run_cli_tests.sh batch tests/rxt/1_sphere nomc

test2_nomc: iad
	IAD_EXECUTABLE=$(IAD_EXECUTABLE) tests/cli/run_cli_tests.sh batch tests/rxt/2_sphere nomc

layertest: $(WSRC) $(NRSRC)
	cd src ; make layer_test
	src/layer_test

wintest: ad.exe iad.exe
	make IAD_EXECUTABLE='wine ./iad.exe' test
	IAD_EXECUTABLE='wine ../iad.exe' tests/cli/run_cli_tests.sh batch

help::
	@echo;\
	echo "Targets available for this Makefile:";\
	echo "  ad            compile forward Adding-Doubling program";\
	echo "  iad           compile inverse Adding-Doubling program";\
	echo "  dist          create a unix distribution";\
	echo "  dists         make unix and windows distributions";\
	echo "  docs          generate source docs, manual, tutorial";\
	echo "  executables   build all the binaries";\
	echo "  install       install ad and iad programs";\
	echo "  install-lib   install interface and library programs";\
	echo "  lib           create library binary and interface files";\
	echo "  windist       create a windows distribution";\
	echo "CLEANING";\
	echo "  clean-generated-dry-run  list generated files that would be removed";\
	echo "  clean-generated          remove ignored build/doc/test outputs";\
	echo "  clean-generated-all      remove generated outputs and tracked generated deliverables";\
	echo "  clean-generated-archives remove generated release archives";\
	echo "  scratch-build            remove generated deliverables, then rebuild tidy/docs/all";\
	echo "TESTING";\
	echo "  longtest      run iad program on a bunch of test files";\
	echo "  shorttest     same as test below";\
	echo "  test0         run tests/rxt/0_sphere files";\
	echo "  test1         run tests/rxt/1_sphere files";\
	echo "  test2         run tests/rxt/2_sphere files";\
	echo "  test1_nomc    run tests/rxt/1_sphere files with -M 0";\
	echo "  test2_nomc    run tests/rxt/2_sphere files with -M 0";\
	echo "  veryshorttest github action test target";\
	echo "  wintest       windows command-line and file tests using wine";\
	echo "MAINTENANCE";\
	echo "  clean         remove most generated objects";\
	echo "  help          this message";\
	echo "  realclean     remove all generated objects";\
	echo "  tidy          generate and format .c and .h files";\


.SECONDARY: $(HSRC) $(CSRC)

.PHONY: clean realclean clean-generated-dry-run clean-generated \
        clean-generated-all-dry-run clean-generated-all scratch-build \
        clean-generated-archives-dry-run clean-generated-archives \
        dists docs test lib install tidy dist windist \
        test veryshorttest shorttest longtest test0 test1 test2 test1_nomc test2_nomc \
        layertest wintest
