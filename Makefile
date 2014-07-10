#
#  Makefile by Scott Prahl, Apr 2013
#

VERSION = 3-9-11

#Base directory for installation
DESTDIR=/usr/local

#typical install hierarchy
BIN_INSTALL=$(DESTDIR)/bin
LIB_INSTALL=$(DESTDIR)/lib
INC_INSTALL=$(DESTDIR)/include

#Library extension
LIB_EXT = .dylib		#MacOS X shared lib
#LIB_EXT = .so			#linux shared lib
#LIB_EXT = .a			#static lib

#Mathematica install requires binaries to be installed and then the
#files stored in the right place
MMA_INSTALL=/Users/prahl/Library/Mathematica/Applications/Optics

TAR=gnutar

CFLAGS = 
CFLAGS = -dynamic -fno-common -Wall -ansi  #Needed flags for MacOS X to build .dylib
LOCAL_LIBRARIES = -lm

IAD_OBJ = src/iad_util.o  src/iad_calc.o src/iad_find.o src/iad_pub.o  src/iad_io.o 

AD_OBJ	= src/nr_zbrak.o  src/ad_bound.o src/ad_doubl.o src/ad_frsnl.o src/ad_globl.o src/ad_matrx.o \
		  src/ad_phase.o  src/ad_prime.o src/ad_radau.o src/ad_start.o src/ad_cone.o  src/ad_layers.o
		  
NR_OBJ	= src/nr_amoeb.o  src/nr_amotr.o src/nr_brent.o src/nr_gaulg.o src/nr_mnbrk.o src/nr_rtsaf.o\
		  src/nr_util.o   src/nr_hj.o

MAIN = Makefile INSTALL README

DOCS =	doc/CHANGELOG        doc/ToDo                  doc/manual.tex      \
		doc/ad_src.pdf       doc/iad_src.pdf	       doc/manual.pdf      \
		doc/ch3RTcorr.pdf    doc/ch3spheremeas.pdf     doc/ch3spheresR.pdf \
		doc/ch3spheresT.pdf  doc/ch3Doublespheres.pdf  doc/colltrans.pdf   \
		doc/lightloss.pdf    doc/niek_graph.pdf        doc/glass-slide.pdf \
		doc/cmdexe.png       doc/valid.png             doc/dual8.png \
        doc/dual90.png       doc/iad.bib

TEST =	\
		test/Makefile       test/basic-A.rxt    test/basic-B.rxt      test/basic-C.rxt    test/basic-D.rxt\
		test/double.rxt     test/example2.rxt   test/il-A.rxt         test/il-B.rxt       test/il-C.rxt\
		test/ink-A.rxt      test/ink-B.rxt      test/ink-C.rxt        test/kenlee-A.rxt   test/kenlee-B.rxt\
		test/kenlee-C.rxt   test/newton.rxt     test/royston2.rxt     test/royston3-A.rxt test/royston3-B.rxt\
		test/royston3-C.rxt test/royston3-D.rxt test/royston3-E.rxt   test/royston9-A.rxt test/royston9-B.rxt\
		test/royston9-C.rxt test/royston9-D.rxt test/royston_ink1.rxt test/sample-A.rxt   test/sample-B.rxt\
		test/sample-C.rxt   test/sample-D.rxt   test/sample-E.rxt     test/sample-F.rxt   test/sample-G.rxt\
		test/sevick-A.rxt   test/sevick-B.rxt   test/terse-A.rxt      test/terse-B.rxt    test/tio2_vis.rxt\
		test/uterus.rxt     test/valid.bat      test/vio-A.rxt        test/vio-B.rxt      test/x_bad_data.rxt\
		test/ville1.rxt     test/fairway-A.rxt  test/fairway-B.rxt    test/fairway-C.rxt  test/fairway-D.rxt\
		test/fairway-E.rxt  test/basic-E.rxt
		
		
MMA =  mma/AD.m mma/AD.tm mma/Makefile

WSRC =	src/ad.w			src/ad_frsnl.w		 src/ad_prime.w	  src/iad_io.w \
		src/ad_globl.w	    src/ad_radau.w	     src/iad_main.w \
		src/ad_bound.w		src/ad_layers.w		 src/ad_start.w	  src/iad_main_mus.w \
		src/ad_chapter.w	src/ad_layers_test.w src/iad.w		  src/iad_pub.w \
		src/ad_cone.w		src/ad_main.w		 src/iad_calc.w	  src/iad_type.w \
		src/ad_cone_test.w	src/ad_matrx.w		 src/iad_util.w	  src/ad_oblique_test.w\
		src/ad_doubl.w		src/ad_phase.w		 src/iad_find.w

NRSRC = src/nr_amoeb.c		src/nr_amotr.h		 src/nr_gaulg.c	  src/nr_mnbrk.h  src/nr_util.c	   src/nr_zbrak.h \
		src/nr_amoeb.h		src/nr_brent.c		 src/nr_gaulg.h	  src/nr_rtsaf.c  src/nr_util.h \
		src/nr_amotr.c		src/nr_brent.h		 src/nr_mnbrk.c	  src/nr_rtsaf.h  src/nr_zbrak.c \
		src/mygetopt.c		src/mygetopt.h		 src/version.h	  src/mc_lost.c	  src/mc_lost.h	  \
		src/nr_hj.c         src/nr_hj.h          src/nr_zbrent.h  src/nr_zbrent.c

CSRC  = src/ad_frsnl.c		src/ad_globl.c		 src/ad_matrx.c	  src/ad_start.c  src/iad_calc.c   src/iad_io.c	  \
		src/iad_main.c		src/ad_doubl.c	     src/iad_util.c   src/ad_radau.c  src/iad_pub.c \
		src/ad_prime.c		src/iad_find.c		 src/ad_phase.c	  src/ad_main.c	  src/ad_bound.c   src/ad_cone.c \
		src/ad_layers.c		src/version.c		 
 
HSRC  = src/ad_bound.h		src/ad_globl.h		 src/ad_phase.h	  src/ad_start.h  src/iad_io.h	   src/iad_type.h \
		src/ad_doubl.h		src/ad_main.h		 src/ad_prime.h	  src/iad_calc.h  src/iad_main.h   src/iad_util.h \
		src/ad_frsnl.h		src/ad_matrx.h		 src/ad_radau.h	  src/iad_find.h  src/iad_pub.h	   src/ad_cone.h  \
		src/ad_cone_ez.h	src/ad_layers.h

OSRC  = src/system.bux src/ad.bux src/iad.bux src/cobweb.pl src/version.pl src/Makefile src/iad.vcproj src/toDOS.pl

#tricky perl script that avoids variable names
#there must be a better way to do this name mangling!
#this is necessary because you can only have one period in a zip archive name

XP_VERSION := $(shell echo $(VERSION) | perl -pe 's/\./-/g; print "xp-version-"')

all :  ad iad
	
check :
	cd test ; make

doc docs: ad_doc iad_doc manual_doc

lib :
	cd src ; make libiad.h
	cp src/libiad.h .
	cd src ; make libiad$(LIB_EXT)
	cp src/libiad$(LIB_EXT) .
	
mma: mma/AD.m mma/AD.exe

mma/AD.exe:
	make install-lib
	cd mma ; make
	
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
	
install-mma: mma/AD.m mma/AD.exe
	make install-lib
	mkdir -p $(MMA_INSTALL)
	mkdir -p $(MMA_INSTALL)/External
	cp mma/AD.m $(MMA_INSTALL)
	cp mma/AD.exe $(MMA_INSTALL)/External
	
install-all:
	make install
	make install-lib
	make install-mma

dist: 
	touch src/version.h
	cd src && ./version.pl
	make doc
	make clean
	make
	cd src ; make tidy
	mkdir -p	   iad-$(VERSION)
	mkdir -p	   iad-$(VERSION)/doc
	mkdir -p	   iad-$(VERSION)/test
	mkdir -p	   iad-$(VERSION)/src
	mkdir -p	   iad-$(VERSION)/mma
	ln $(MAIN)	   iad-$(VERSION)
	ln $(DOCS)	   iad-$(VERSION)/doc
	ln $(TEST)	   iad-$(VERSION)/test
	ln $(WSRC)	   iad-$(VERSION)/src
	ln $(HSRC)	   iad-$(VERSION)/src
	ln $(CSRC)	   iad-$(VERSION)/src
	ln $(NRSRC)	   iad-$(VERSION)/src
	ln $(OSRC)	   iad-$(VERSION)/src
	ln $(MMA)	   iad-$(VERSION)/mma
#	$(TAR) cvf - iad-$(VERSION) | gzip > iad-$(VERSION).tar.gz
	zip -r iad-$(VERSION) iad-$(VERSION)
	rm -rf iad-$(VERSION)

zip: 
#	make doc
#	cd src ; make tidy
	mkdir -p	       $(XP_VERSION)
	mkdir -p	       $(XP_VERSION)/doc
	mkdir -p	       $(XP_VERSION)/test
	mkdir -p	       $(XP_VERSION)/src
	ln INSTALL	       $(XP_VERSION)
	ln README	       $(XP_VERSION)
	ln iad.exe	       $(XP_VERSION)
	ln ad.exe	       $(XP_VERSION)
	ln doc/manual.pdf  $(XP_VERSION)/doc
	ln doc/iad_src.pdf $(XP_VERSION)/doc
	ln doc/ad_src.pdf  $(XP_VERSION)/doc
	ln $(TEST)	       $(XP_VERSION)/test
	ln $(HSRC)	       $(XP_VERSION)/src
	ln $(CSRC)	       $(XP_VERSION)/src
	ln $(NRSRC)	       $(XP_VERSION)/src
	ln src/iad.vcproj  $(XP_VERSION)/src
	rm $(XP_VERSION)/src/ad_cone.c    
	rm $(XP_VERSION)/src/ad_cone.h 
	rm $(XP_VERSION)/src/ad_cone_ez.h
	rm $(XP_VERSION)/src/ad_layers.c  
	rm $(XP_VERSION)/src/ad_layers.h
	rm $(XP_VERSION)/src/ad_main.c
	rm $(XP_VERSION)/src/ad_main.h
	rm $(XP_VERSION)/src/nr_hj.c      
	rm $(XP_VERSION)/src/nr_hj.h
	src/toDOS.pl $(XP_VERSION)/src/*.c
	src/toDOS.pl xp_iad/src/*.h
	src/toDOS.pl $(XP_VERSION)/test/*.rxt
	src/toDOS.pl $(XP_VERSION)/test/valid.bat
	rm $(XP_VERSION)/src/*.bak 
	rm $(XP_VERSION)/test/*.bak
	rm $(XP_VERSION)/test/Makefile
	zip -r $(XP_VERSION) $(XP_VERSION)
	rm -rf $(XP_VERSION)

ad: $(WSRC) $(NRSRC)
	cd src ; make ad
	cp src/ad ad$(BINARY_EXTENSION)

iad: $(WSRC) $(NRSRC)
	cd src ; make iad
	cp src/iad iad$(BINARY_EXTENSION)

ad_debug: $(WSRC) $(NRSRC)
	cd src ; make ad_debug
	cp src/ad ad$(BINARY_EXTENSION)

iad_debug: $(WSRC) $(NRSRC)
	cd src ; make iad_debug
	cp src/iad iad$(BINARY_EXTENSION)

ad_doc: 
	perl -pi -e 's/\\def\\adversion.*/\\def\\adversion{$(VERSION)}/' src/ad.w
	cd src ; make ad_doc

iad_doc:
	perl -pi -e 's/\\def\\iadversion.*/\\def\\iadversion{$(VERSION)}/' src/iad.w
	cd src ; make iad_doc

manual_doc doc/manual.pdf: doc/manual.tex 
	cd doc ; pdflatex manual
	cd doc ; bibtex manual
	cd doc ; pdflatex manual
	cd doc ; pdflatex manual
	
clean:
	rm -f ad iad *.pdf libiad.a libiad.so libiad.dylib libiad.h src/lib_iad.h src/lib_ad.h
	rm -f src/*.o test/*.abg
	rm -f texexec.tmp texexec.tui mpgraph.mp mprun.mp texexec-mpgraph.mp texexec.tex
	rm -f src/*.aux src/*.dvi src/*.idx src/*.ref src/*.sref src/*.tex src/*.toc src/*.log src/*.scn
	
realclean:
	make clean
	cd src ; make realclean
	rm -f mpgraph.mp texexec-mpgraph.mp texexec.tmp
	rm -f ad iad libiad.h libiad$(LIB_EXT)
	rm -f doc/manual.pdf doc/manual.bbl doc/manual.blg doc/manual.out doc/manual.toc 
	cd mma	; make clean
	cd test ; make clean
	rm -f $(CSRC) $(HSRC)
	
test:
	@echo "********* Basic tests ***********"
	./iad -V 0 -r 0
	./iad -V 0 -r 1
	./iad -V 0 -r 0.4
	./iad -V 0 -r 0.4 -t 0.1
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002
	./iad -V 0 -r 0.4 -t 0.1 -u 0.049787
	@echo "********* Specify sample index ************"
	./iad -V 0 -r 0.4 -n 1.5
	./iad -V 0 -r 0.4 -t 0.1 -n 1.5
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.5
	./iad -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.5
	@echo "********* Specify slide index ************"
	./iad -V 0 -r 0.4 -n 1.4 -N 1.5
	./iad -V 0 -r 0.4 -t 0.1 -n 1.4 -N 1.5
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.4 -N 1.5
	./iad -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5
	@echo "********* One slide on top ************"
	./iad -V 0 -r 0.4 -n 1.4 -N 1.5 -G t
	./iad -V 0 -r 0.4 -t 0.1 -n 1.4 -N 1.5 -G t
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.4 -N 1.5 -G t
	./iad -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5 -G t
	@echo "********* One slide on bottom ************"
	./iad -V 0 -r 0.4 -n 1.4 -N 1.5 -G b
	./iad -V 0 -r 0.4 -t 0.1 -n 1.4 -N 1.5 -G b
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.4 -N 1.5 -G b
	./iad -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5 -G b
	@echo "********* Absorbing Slide Tests ***********"
	./iad -V 0 -r 0.0000000 -t 0.135335 -E 0.5
	./iad -V 0 -r 0.0249268 -t 0.155858 -E 0.5
	./iad -V 0 -r 0.0399334 -t 0.124727 -E 0.5 -n 1.5 -N 1.5
	./iad -V 0 -r 0.0520462 -t 0.134587 -E 0.5 -n 1.5 -N 1.5
	@echo "********* Constrain g ************"
	./iad -V 0 -r 0.4        -g 0.9
	./iad -V 0 -r 0.4 -t 0.1 -g 0.9
	./iad -V 0 -r 0.4        -g 0.9 -n 1.5
	./iad -V 0 -r 0.4 -t 0.1 -g 0.9 -n 1.5
	./iad -V 0 -r 0.4        -g 0.9 -n 1.4 -N 1.5
	./iad -V 0 -r 0.4 -t 0.1 -g 0.9 -n 1.4 -N 1.5
	@echo "********* Constrain a ************"
	./iad -V 0 -r 0.4        -a 0.9
	./iad -V 0 -r 0.4 -t 0.1 -a 0.9
	@echo "********* Constrain b ************"
	./iad -V 0 -r 0.4        -b 3
	./iad -V 0 -r 0.4 -t 0.1 -b 3
	./iad -V 0 -r 0.4 -t 0.1 -b 3 -n 1.5
	./iad -V 0 -r 0.4 -t 0.1 -b 3 -n 1.4 -N 1.5
	@echo "********* Constrain mu_s ************"
	./iad -V 0 -r 0.4        -F 30
	./iad -V 0 -r 0.4 -t 0.1 -F 30
	./iad -V 0 -r 0.4        -F 30 -n 1.5
	./iad -V 0 -r 0.4 -t 0.1 -F 30 -n 1.5
	./iad -V 0 -r 0.4        -F 30 -n 1.4 -N 1.5
	./iad -V 0 -r 0.4 -t 0.1 -F 30 -n 1.4 -N 1.5
	@echo "********* Constrain mu_a ************"
	./iad -V 0 -r 0.3        -A 0.6
	./iad -V 0 -r 0.3 -t 0.1 -A 0.6
	./iad -V 0 -r 0.3        -A 0.6 -n 1.5
	./iad -V 0 -r 0.3 -t 0.1 -A 0.6 -n 1.5
	./iad -V 0 -r 0.3        -A 0.6 -n 1.4 -N 1.5
	./iad -V 0 -r 0.3 -t 0.1 -A 0.6 -n 1.4 -N 1.5
	@echo "********* Constrain mu_a and g************"
	./iad -V 0 -r 0.3        -A 0.6 -g 0.6
	./iad -V 0 -r 0.3        -A 0.6 -g 0.6 -n 1.5
	./iad -V 0 -r 0.3        -A 0.6 -g 0.6 -n 1.4 -N 1.5
	@echo "********* Constrain mu_s and g************"
	./iad -V 0 -r 0.3        -F 2.0 -g 0.5
	./iad -V 0 -r 0.3        -F 2.0 -g 0.5 -n 1.5
	./iad -V 0 -r 0.3        -F 2.0 -g 0.5 -n 1.4 -N 1.5
	@echo "********* Basic One Sphere tests ***********"
	./iad -V 0 -r 0.4                     -S 1
	./iad -V 0 -r 0.4 -t 0.1              -S 1
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002     -S 1
	./iad -V 0 -r 0.2 -t 0.2 -u 0.0049787 -S 1
	@echo "******** Basic 10,000 photon tests *********"
	./iad -V 0 -r 0.4                     -S 1 -p 10000
	./iad -V 0 -r 0.4 -t 0.1              -S 1 -p 10000
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002     -S 1 -p 10000
	./iad -V 0 -r 0.2 -t 0.2 -u 0.0049787 -S 1 -p 10000
	@echo "******** Basic timed photon tests *********"
	./iad -V 0 -r 0.4                     -S 1 -p -500
	./iad -V 0 -r 0.4 -t 0.1              -S 1 -p -500
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002     -S 1 -p -500
	./iad -V 0 -r 0.2 -t 0.2 -u 0.0049787 -S 1 -p -500
	@echo "********* More One Sphere tests ***********"
	./iad -V 0 -r 0.4                     -S 1 -1 '200 13 13 0'
	./iad -V 0 -r 0.4 -t 0.1              -S 1 -1 '200 13 13 0'
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002     -S 1 -1 '200 13 13 0'
	./iad -V 0 -r 0.2 -t 0.2 -u 0.0049787 -S 1 -1 '200 13 13 0'
	@echo "********* Basic Two Sphere tests ***********"
	./iad -V 0 -r 0.4                     -S 2
	./iad -V 0 -r 0.4 -t 0.1              -S 2
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002     -S 2
	./iad -V 0 -r 0.2 -t 0.1 -u 0.0049787 -S 2
	@echo "********* More Two Sphere tests ***********"
	./iad -V 0 -r 0.4                     -S 2 -1 '200 13 13 0' -2 '200 13 13 0'
	./iad -V 0 -r 0.4 -t 0.1              -S 2 -1 '200 13 13 0' -2 '200 13 13 0'
	./iad -V 0 -r 0.4 -t 0.1 -u 0.002     -S 2 -1 '200 13 13 0' -2 '200 13 13 0'
	./iad -V 0 -r 0.2 -t 0.1 -u 0.0049787 -S 2 -1 '200 13 13 0' -2 '200 13 13 0'
	@echo "********* Oblique tests ***********"
	./iad -V 0 -i 60 -r 0.00000 -t 0.13691
	./iad -V 0 -i 60 -r 0.14932 -t 0.23181
	./iad -V 0 -i 60 -r 0.61996 -t 0.30605
	@echo "********* Different Tstd and Rstd ***********"
	./iad -V 0 -r 0.4 -t 0.005  -d 1 -M 0  -S 1 -1 '200 13 13 0' -T 0.5
	./iad -V 0 -r 0.4 -t 0.005  -d 1 -M 0  -S 1 -1 '200 13 13 0'
	./iad -V 0 -r 0.2 -t 0.01   -d 1 -M 0  -S 1 -1 '200 13 13 0' -R 0.5
	./iad -V 0 -r 0.2 -t 0.01   -d 1 -M 0  -S 1 -1 '200 13 13 0'
	
help::
	@echo;\
	echo "Targets available for this Makefile:";\
	echo "ad		  -- compile forward Adding-Doubling program";\
	echo "iad		  -- compile inverse Adding-Doubling program";\
	echo "clean		  -- remove most generated objects";\
	echo "docs		  -- generate TEX out of all files";\
	echo "dist		  -- gzipped tarball of all files";\
	echo "install	  -- install ad and iad programs";\
	echo "install-lib -- install interface and library programs";\
	echo "install-mma -- install Mathematica files";\
	echo "install-all -- install all the above";\
	echo "lib		  -- library binary and interface files";\
	echo "mma		  -- adding-doubling files for Mathematica interface";\
	echo "realclean	  -- remove all generated objects";\
	echo "test		  -- run iad program on a bunch of test files";\
	echo "tidy		  -- generate .c and .h files"

.SECONDARY: $(HSRC) $(CSRC)

.PHONY : clean realclean dist docs test check lib mma install install-all iad_doc ad_doc manual_doc new_version
