all : ad iad

INDENT = indent

INDENT_OPT = -bad -bap -c33 -ci1 -d1 -di20 -i4 -l78 -nbc -nfc1 -psl
INDENT_OPT =

ASSERT = -DNDEBUG
#ASSERT =

LOCAL_LIBRARIES = -lm
SHARED_LIB_OPT = -lm -shared 
DYNAMIC_LIB_OPT = -dynamiclib -install_name /usr/local/lib/libiad.dylib \
				  -compatibility_version 2.0 -current_version 2.0.0 \

CFLAGS = -fno-common -Wall $(ASSERT) -pedantic -g -fPIC -fsigned-char

IAD_OBJ = iad_util.o  iad_calc.o iad_find.o iad_pub.o  iad_io.o 

AD_OBJ	= nr_zbrak.o  ad_bound.o ad_doubl.o ad_frsnl.o ad_globl.o ad_matrx.o \
		  ad_phase.o  ad_prime.o ad_radau.o ad_start.o ad_cone.o  ad_layers.o
		  
NR_OBJ	= nr_amoeb.o  nr_amotr.o nr_brent.o nr_gaulg.o nr_mnbrk.o nr_rtsaf.o\
		  nr_util.o	  mc_lost.o  nr_hj.o    nr_zbrent.o

ad: $(AD_OBJ) $(NR_OBJ) ad_main.o mygetopt.o version.o
	./version.pl 
	$(CC) -c -o version.o version.c 	
	$(CC) -o $@ $^ $(LOCAL_LIBRARIES)

iad: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) iad_main.o mygetopt.o version.o
	./version.pl 
	$(CC) -c -o version.o version.c 
	$(CC) -o $@ $^ $(LOCAL_LIBRARIES)

libiad.a : $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) libiad.h version.o
	$(AR) cr libiad.a $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o
	ranlib libiad.a

libiad.so: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) libiad.h version.o
	$(CC) $(SHARED_LIB_OPT) -o libiad.so $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o

libiad.dylib: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) libiad.h version.o
	$(CC) $(DYNAMIC_LIB_OPT) -o libiad.dylib $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o

libiad.dll: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) libiad.h version.o
	$(CC) -shared -o libiad.dll $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o

libiad.h : ad_prime.w iad_pub.w ad_cone.w ad_layers.w
	cd src ; ctangle -bhp ad_prime
	cd src ; ctangle -bhp iad_pub
	cd src ; ctangle -bhp ad_cone
	cd src ; ctangle -bhp ad_layers
	cat lib_ad.h lib_iad.h ad_cone_ez.h ad_layers.h > libiad.tmp
	./cobweb.pl < libiad.tmp > libiad.tmp2
	$(INDENT) $(INDENT_OPT) < libiad.tmp2 > libiad.h
	rm -f libiad.tmp lib_ad.h lib_iad.h libiad.tmp2

ad_cone_test: $(AD_OBJ) $(NR_OBJ) ad_cone_test.o mygetopt.o version.o
	./version.pl 
	$(CC) -c -o version.o version.c 	
	$(CC) -o $@ $^ $(LOCAL_LIBRARIES)

ad_oblique_test: $(AD_OBJ) $(NR_OBJ) ad_oblique_test.o mygetopt.o version.o
	./version.pl 
	$(CC) -c -o version.o version.c 	
	$(CC) -o $@ $^ $(LOCAL_LIBRARIES)

ad_doc: ad.w
	cweave -bph ad || (echo "cweave exited with code $$? - continue anyway"; exit 0)
	pdftex ad
	qpdf ad.pdf --pages ad.pdf 129-130 ad.pdf 1-128 -- ad2.pdf
	mv ad2.pdf ../doc/ad_src.pdf
	rm ad.pdf 

#	cd src ; ctwill -bhp ad.w
#	cd src ; ctwill -bhp ad.w
#	tex -output-directory src ad.tex
#	cd src ; refsort < ad.ref > ad.sref
#	pdftex -output-directory src src/ad.tex
#	texexec --pdfselect --paper=letter --selection=1:132 --result src/ad_body.pdf src/ad.pdf
#	cd src ; cweave -bhp ad
#	pdftex -output-directory src src/ad.tex
#	texexec --pdfselect --paper=letter --selection=113:123 --result src/ad_index.pdf src/ad.pdf
#	texexec --pdfselect --paper=letter --selection=124:125 --result src/ad_toc.pdf	 src/ad.pdf
#	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=doc/ad_src.pdf src/ad_toc.pdf src/ad_body.pdf src/ad_index.pdf &> /dev/null 
#	rm -f src/ad_body.pdf src/ad_toc.pdf src/ad_index.pdf src/ad.pdf

iad_doc: iad.w
	cweave -bph iad || (echo "cweave exited with code $$? - continue anyway"; exit 0)
	pdftex iad
	qpdf iad.pdf --pages iad.pdf 142 iad.pdf 1-141 -- iad2.pdf
	mv iad2.pdf ../doc/iad_src.pdf
	rm iad.pdf 
#	cd src ; ctwill -bhp iad.w
#	cd src ; ctwill -bhp iad.w
#	cd src ; tex iad.tex
#	cd src ; refsort < iad.ref > iad.sref
#	pdftex -output-directory src src/iad.tex
#	texexec --pdfselect --paper=letter --selection=1:113 --result src/iad_body.pdf src/iad.pdf
#	cd src ; cweave -bhp iad
#	pdftex	-output-directory src src/iad.tex
#	texexec --pdfselect --paper=letter --selection=96:100 --result src/iad_index.pdf src/iad.pdf
#	texexec --pdfselect --paper=letter --selection=106	  --result src/iad_toc.pdf	 src/iad.pdf
#	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=doc/iad_src.pdf src/iad_toc.pdf src/iad_body.pdf src/iad_index.pdf &> /dev/null 
#	rm -f src/iad_body.pdf src/iad_toc.pdf src/iad_index.pdf src/iad.pdf
	
CSRC  = src/ad_frsnl.c		src/ad_globl.c		 src/ad_matrx.c	  src/ad_start.c  src/iad_calc.c   src/iad_io.c	  \
		src/iad_type.c		src/iad_main.c		 src/ad_doubl.c	  src/iad_util.c  src/ad_radau.c   src/iad_pub.c \
		src/ad_prime.c		src/iad_find.c		 src/ad_phase.c	  src/ad_main.c	  src/ad_bound.c   src/ad_cone.c \
		src/ad_layers.c		src/version.c		 

tidy:
	apply ctangle *.w
	./version.pl
	apply ./cobweb.pl ad*.c ad*.h iad*.c iad*.h
	rm -f iad.c iad_main_mus.c
	rm -f ad_boil.c
	rm -f *.bak *.BAK
	rm -f iad_type.c
	rm -f iad_main.h
	gindent $(INDENT_OPT) ad*.c iad*.c ad*.h iad*.h
	rm *~

ad_oblique_test: $(AD_OBJ) $(NR_OBJ) ad_oblique_test.o mygetopt.o version.o

layer_test: $(AD_OBJ) $(NR_OBJ) mygetopt.o version.o ad_layers.o ad_layers_test.o	
	ctangle -bhp ad_layers.w
	$(CC) -c -o ad_layers.o ad_layers.c 
	ctangle -bhp ad_layers_test.w
	$(CC) -c -o ad_layers_test.o ad_layers_test.c 
	./version.pl 
	$(CC) -c -o version.o version.c 	
	$(CC) -o $@ $^ $(LOCAL_LIBRARIES)

# non-trivial dependencies
lib_ad.h   :  ad
lib_iad.h  :  iad

ad_globl.o :  ad_frsnl.o	
ad_bound.o :  ad_globl.o ad_matrx.o ad_frsnl.o	 
ad_doubl.o :  ad_globl.o nr_util.o	ad_matrx.o
ad_start.o :  ad_globl.o nr_util.o	ad_frsnl.o nr_gaulg.o ad_radau.o ad_phase.o \
				  ad_frsnl.o ad_matrx.o ad_phase.o
ad_cone.o  :  ad_cone.w ad_globl.o nr_util.o ad_bound.o ad_doubl.o ad_start.o
ad_prime.o :  ad_prime.w ad_cone.o
ad_main.o  :  ad_main.w ad_globl.o ad_prime.o mygetopt.o mc_lost.o

iad_calc.o :  iad_calc.w ad_globl.o ad_frsnl.o ad_prime.o iad_type.h iad_util.o \
              nr_zbrent.o
iad_pub.o  :  iad_pub.w ad_globl.o nr_util.o iad_type.h iad_find.o iad_io.o
iad_util.o :  iad_util.w ad_globl.o nr_util.o ad_frsnl.o ad_bound.o iad_type.h iad_calc.h iad_pub.h
iad_find.o :  iad_find.w iad_find.h ad_globl.o nr_util.o iad_type.h iad_util.o iad_calc.o
iad_main.o :  iad_main.w ad_globl.o iad_type.h iad_pub.o iad_io.o mygetopt.o iad_calc.o
iad_io.o   :  iad_io.w iad_type.h iad_pub.h 
nr_amoeb.o :  iad_type.h 
ad_layers.o: nr_util.o ad_globl.o ad_bound.o ad_doubl.o ad_prime.o ad_matrx.o ad_prime.o

ad_bound.c ad_bound.h : ad_bound.w
	ctangle -bhp $<

ad_cone.c ad_cone.h: ad_cone.w
	ctangle -bhp $<

ad_doubl.c ad_doubl.h: ad_doubl.w
	ctangle -bhp $<

ad_frsnl.c ad_frsnl.h: ad_frsnl.w
	ctangle -bhp $<

ad_globl.c ad_globl.h: ad_globl.w
	ctangle -bhp $<

ad_layers.c ad_layers.h: ad_layers.w
	ctangle -bhp $<

ad_main.c ad_main.h: ad_main.w
	ctangle -bhp $<

ad_matrx.c ad_matrx.h: ad_matrx.w
	ctangle -bhp $<

ad_prime.c ad_prime.h: ad_prime.w
	ctangle -bhp $<

ad_radau.c ad_radau.h: ad_radau.w
	ctangle -bhp $<

ad_start.c ad_start.h: ad_start.w
	ctangle -bhp $<

iad_calc.c iad_calc.h: iad_calc.w
	ctangle -bhp $<

iad_find.c iad_find.h: iad_find.w
	ctangle -bhp $<

iad_io.c iad_io.h: iad_io.w
	ctangle -bhp $<

iad_main.c iad_main.h: iad_main.w
	ctangle -bhp $<

iad_pub.c iad_pub.h: iad_pub.w
	ctangle -bhp $<

iad_type.h : iad_type.w
	ctangle -bhp $<
	rm iad_type.c

iad_util.c iad_util.h: iad_util.w
	ctangle -bhp $<

mc_lost.c : iad_type.h

version.o : 
	./version.pl 
	$(CC) -c -o version.o version.c 

clean:
	rm -f $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) iad_main.o 
	rm -f mygetopt.o version.o ad_main.o ad_cone_test.o ad_oblique_test.o
	
realclean:
	make clean
	rm -f ad*.h ad*.c iad*.h iad*.c libiad.h
	rm -f libiad.* lib_ad* lib_iad*
	rm -f ad.idx ad.log ad.scn ad.tex ad.toc ad.pdf
	rm -f iad.idx iad.log iad.scn iad.tex iad.toc iad.pdf
	
.PHONEY : clean realclean lib

.PRECIOUS: %.c %.h
