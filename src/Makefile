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

CFLAGS = -fno-common -Wall $(ASSERT) -pedantic -g -fPIC

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
	cweave -bph ad
	pdftex ad
	qpdf ad.pdf --pages ad.pdf 121-122 ad.pdf 1-120 -- ad2.pdf
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
	cweave -bph iad
	pdftex iad
#	mv iad.pdf ../doc/iad_src.pdf
	qpdf iad.pdf --pages iad.pdf 131 iad.pdf 1-130 -- iad2.pdf
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
iad_find.o :  iad_find.w ad_globl.o nr_util.o iad_type.h iad_util.o iad_calc.o
iad_main.o :  iad_main.w ad_globl.o iad_type.h iad_pub.o iad_io.o mygetopt.o iad_calc.o
iad_io.o   :  iad_io.w iad_type.h iad_pub.h 
nr_amoeb.o :  iad_type.h 

iad_calc.h :  iad_calc.w		#needed because iad_util.o includes iad_calc.h
	ctangle -bhp iad_calc.w

iad_type.h : iad_type.w
	ctangle -bhp iad_type.w
	rm iad_type.c

iad_pub.h : iad_pub.w          #needed because iad_io.o includes iad_pub.h
	ctangle -bhp iad_pub.w

ad_bound.h : ad_bound.w
	ctangle -bhp ad_bound.w
	
mc_lost.c : iad_type.h

version.o : 
	./version.pl 
	$(CC) -c -o version.o version.c 

clean:
	rm -f $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) iad_main.o 
	rm -f mygetopt.o version.o ad_main.o ad_cone_test.o ad_oblique_test.o
	
realclean:
	make clean
	rm -f ad*.h ad*.c iad*.h iad*.c
	rm -f libiad.* lib_ad* lib_iad*
	rm -f ad.idx ad.log ad.scn ad.tex ad.toc ad.pdf
	rm -f iad.idx iad.log iad.scn iad.tex iad.toc iad.pdf
	
.PHONEY : clean realclean lib
