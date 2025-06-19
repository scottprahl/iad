all : ad iad

INDENT = indent
INDENT_OPT = -br -brs -nbs -ci4 -cli0 -bli0 -l120 -ncdb -nce -nlp -npcs -npsl -i4 -nut

ASSERT =
ASSERT = -DNDEBUG

SHARED_LIB_OPT = -lm -shared
DYNAMIC_LIB_OPT = -dynamiclib -install_name /usr/local/lib/libiad.dylib \
                  -compatibility_version 2.0 -current_version 2.0.0 \

CFLAGS = -fno-common -Wall $(ASSERT) -pedantic -g -fPIC -fsigned-char

IAD_OBJ = iad_util.o  iad_calc.o iad_find.o iad_pub.o  iad_io.o   cache.o  mc_lost.o

AD_OBJ  = nr_zbrak.o  ad_bound.o ad_doubl.o ad_frsnl.o ad_globl.o ad_matrx.o \
          ad_phase.o  ad_prime.o ad_radau.o ad_start.o ad_cone.o  ad_layers.o

NR_OBJ  = nr_amoeb.o  nr_amotr.o nr_brent.o nr_gaulg.o  nr_mnbrk.o     nr_rtsaf.o \
          nr_util.o   nr_hj.o    nr_zbrent.o

cache.o: cache.c cache.h
	$(CC) -c -o cache.o cache.c

ad: $(AD_OBJ) $(NR_OBJ) ad_main.o version.o
	./version.pl
	$(CC) -c -o version.o version.c
	$(CC) -o $@ $^ -lm

iad: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) iad_main.o version.o
	./version.pl
	$(CC) -c -o version.o version.c
	$(CC) -o $@ $^ -lm

libiad.a : $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) libiad.h version.o
	$(AR) cr libiad.a $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o
	ranlib libiad.a

libiad.so: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) libiad.h version.o
	$(CC) $(SHARED_LIB_OPT) -o libiad.so $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o

libiad.dylib: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) libiad.h version.o
	$(CC) $(DYNAMIC_LIB_OPT) -o libiad.dylib $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o

libiad.dll: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) libiad.h version.o
	$(CC) -shared -o libiad.dll $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o

libiad.h : ad iad ad_cone_ez.h ad_layers.h
	cat lib_ad.h lib_iad.h ad_cone_ez.h ad_layers.h > libiad.tmp
	./cobweb.pl < libiad.tmp > libiad.tmp2
	$(INDENT) $(INDENT_OPT) < libiad.tmp2 > libiad.h
	rm -f libiad.tmp lib_ad.h lib_iad.h libiad.tmp2

ad_doc: ad.w
	cweave -bph ad || (echo "cweave exited with code $$? - continue anyway"; exit 0)
	pdftex ad
	qpdf ad.pdf --pages ad.pdf r2-r1 ad.pdf 1-r3 -- ad2.pdf
	mv ad2.pdf ../docs/ad_src.pdf
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
#	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=docs/ad_src.pdf src/ad_toc.pdf src/ad_body.pdf src/ad_index.pdf &> /dev/null
#	rm -f src/ad_body.pdf src/ad_toc.pdf src/ad_index.pdf src/ad.pdf

iad_doc: iad.w
	cweave -bph iad || (echo "cweave exited with code $$? - continue anyway"; exit 0)
	pdftex iad
	qpdf iad.pdf --pages iad.pdf r2-r1 iad.pdf 1-r3 -- iad2.pdf
	mv iad2.pdf ../docs/iad_src.pdf
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
#	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=docs/iad_src.pdf src/iad_toc.pdf src/iad_body.pdf src/iad_index.pdf &> /dev/null
#	rm -f src/iad_body.pdf src/iad_toc.pdf src/iad_index.pdf src/iad.pdf

adchapter: $(AD_OBJ) $(NR_OBJ) ad_chapter.o version.o
	./version.pl
	$(CC) -c -o version.o version.c
	$(CC) -o $@ $^ -lm

layer_test: $(AD_OBJ) $(NR_OBJ) version.o ad_layers.o ad_layers_test.o
	ctangle -bhp ad_layers.w
	$(CC) -c -o ad_layers.o ad_layers.c
	ctangle -bhp ad_layers_test.w
	$(CC) -c -o ad_layers_test.o ad_layers_test.c
	./version.pl
	$(CC) -c -o version.o version.c
	$(CC) -o $@ $^ -lm

cone_test: $(AD_OBJ) $(NR_OBJ) ad_cone_test.o version.o
	./version.pl
	$(CC) -c -o version.o version.c
	$(CC) -o $@ $^ -lm

oblique_test: $(AD_OBJ) $(NR_OBJ) ad_oblique_test.o version.o
	./version.pl
	$(CC) -c -o version.o version.c
	$(CC) -o $@ $^ -lm

mc_test: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o mc_lost.o mc_test.o
	$(CC) -o $@ $^ -lm

mc_lost: $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) version.o mc_lost.o mc_lost_main.o
	$(CC) -o $@ $^ -lm

# ctangle rule to make .c and .h files
%.c %.h: %.w
	ctangle -bhp $<

# non-trivial dependencies
ad_globl.o : ad_globl.c ad_globl.h ad_frsnl.h
ad_bound.o : ad_bound.c ad_bound.h ad_globl.h ad_matrx.h ad_frsnl.h nr_util.h
ad_doubl.o : ad_doubl.c ad_doubl.h ad_globl.h ad_matrx.h nr_util.h
ad_start.o : ad_start.c ad_start.h ad_globl.h ad_frsnl.h ad_radau.h ad_phase.h \
             ad_matrx.h nr_util.h nr_gaulg.h
ad_cone.o  : ad_cone.c ad_cone.h ad_globl.h ad_bound.h ad_matrx.h ad_doubl.h \
             ad_start.h nr_util.h
ad_prime.o : ad_prime.c ad_prime.h ad_globl.h ad_bound.h ad_start.h ad_doubl.h \
             ad_matrx.h ad_cone.h
ad_main.o  : ad_main.c ad_globl.h ad_prime.h ad_cone.h version.h
ad_layers.o: ad_layers.c ad_layers.h nr_util.h ad_globl.h ad_bound.h ad_doubl.h \
             ad_prime.h ad_matrx.h

iad_calc.o : iad_calc.c iad_calc.h ad_globl.h ad_frsnl.h ad_prime.h iad_type.h \
             iad_util.h nr_util.h nr_zbrent.h mc_lost.h
iad_pub.o  : iad_pub.c iad_pub.h ad_globl.h ad_frsnl.h iad_type.h iad_util.h iad_calc.h \
             iad_find.h iad_io.h nr_util.h mc_lost.h
iad_util.o : iad_util.c iad_util.h ad_globl.h ad_frsnl.h ad_bound.h iad_type.h \
             iad_calc.h iad_pub.h nr_util.h
iad_find.o : iad_find.c iad_find.h ad_globl.h iad_type.h iad_util.h iad_calc.h \
             nr_mnbrk.h nr_brent.h nr_amoeb.h nr_util.h
iad_main.o : iad_main.c ad_globl.h ad_frsnl.h ad_prime.h iad_type.h iad_pub.h iad_io.h \
             iad_calc.h iad_util.h mc_lost.h
iad_io.o   : iad_io.c iad_io.h ad_globl.h iad_type.h iad_pub.h version.h

nr_amoeb.o :  iad_type.h
mc_lost.o  :  iad_type.h iad_util.h
cache.o    :  cache.c cache.h

mc_lost_main.o : iad_type.h ad_globl.h ad_prime.h iad_type.h mc_lost.h ad_cone.h
mc_test.o : iad_type.h ad_globl.h ad_prime.h iad_type.h mc_lost.h ad_cone.h

version.o :
	./version.pl
	$(CC) -c -o version.o version.c

clean:
	rm -f $(AD_OBJ) $(IAD_OBJ) $(NR_OBJ) iad_main.o
	rm -f version.o ad_main.o ad_cone_test.o ad_oblique_test.o
	rm -f mc_lost_main.o ad_layers_test.o
	rm -f cone_test layer_test mc_test oblique_test mc_lost
	rm -f ad iad adchapter

realclean:
	make clean
	rm -f ad*.h ad*.c iad*.h iad*.c libiad.h
	rm -f libiad.* lib_ad* lib_iad*
	rm -f ad.idx ad.log ad.scn ad.tex ad.toc ad.pdf
	rm -f iad.idx iad.log iad.scn iad.tex iad.toc iad.pdf

# 'make tidy' does not work when called from the src directory
tidy:
	apply ctangle *.w
	./version.pl
	rm -f iad_main.h ad_main.h iad_type.c iad.c iad_main_mus.c ad_boil.c ad.c
	apply ./cobweb.pl ad*.c ad*.h iad*.c iad*.h lib_ad.h lib_iad.h
	rm -f *.bak *.BAK
	@for file in ad*.c ad*.h iad*.c iad*.h lib_*.h; do \
		{ echo '/* Autogenerated v$(VERSION) from https://github.com/scottprahl/iad */\n'; cat "$$file"; } > temp && mv temp "$$file"; \
	done
	gindent $(INDENT_OPT) ad*.c iad*.c ad*.h iad*.h
	rm *~

.PHONEY : clean realclean lib tidy

.PRECIOUS: %.c %.h
