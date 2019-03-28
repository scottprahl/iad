#
#  Makefile by Scott Prahl in Nov 1998
#

IAD_EXECUTABLE ?= ../iad

TST_FILES = \
basic-A.rxt       fairway-B.rxt     ink-B.rxt         royston3-B.rxt    royston_ink1.rxt  sevick-A.rxt      vio-B.rxt\
basic-B.rxt       fairway-C.rxt     ink-C.rxt         royston3-C.rxt    sample-A.rxt      sevick-b.rxt      x_bad_data.rxt\
basic-C.rxt       fairway-D.rxt     kenlee-A.rxt      royston3-D.rxt    sample-B.rxt      terse-A.rxt\
basic-D.rxt       fairway-E.rxt     kenlee-B.rxt      royston3-E.rxt    sample-C.rxt      terse-B.rxt\
basic-E.rxt       il-A.rxt          kenlee-C.rxt      royston9-A.rxt    sample-D.rxt      tio2_vis.rxt\
double.rxt        il-B.rxt          newton.rxt        royston9-B.rxt    sample-E.rxt      uterus.rxt\
example2.rxt      il-C.rxt          royston2.rxt      royston9-C.rxt    sample-F.rxt      ville1.rxt\
fairway-A.rxt     ink-A.rxt         royston3-A.rxt    royston9-D.rxt    sample-G.rxt      vio-A.rxt

OUT_FILES =  \
basic-A.txt       fairway-B.txt     ink-B.txt         royston3-B.txt    royston_ink1.txt  vio-B.txt\
basic-B.txt       fairway-C.txt     ink-C.txt         royston3-C.txt    sample-A.txt      x_bad_data.txt\
basic-C.txt       fairway-D.txt     kenlee-A.txt      royston3-D.txt    sample-B.txt      terse-A.txt\
basic-D.txt       fairway-E.txt     kenlee-B.txt      royston3-E.txt    sample-C.txt      terse-B.txt\
il-A.txt          kenlee-C.txt      royston9-A.txt    sample-D.txt      tio2_vis.txt \
double.txt        il-B.txt          newton.txt        royston9-B.txt    sample-E.txt      \
example2.txt      il-C.txt          royston2.txt      royston9-C.txt    sample-F.txt      \
fairway-A.txt     ink-A.txt         royston3-A.txt    royston9-D.txt    sample-G.txt      vio-A.txt\

SPECIAL_OUT_FILES = ville1.txt basic-E.txt sevick-A.txt sevick-b.txt  uterus.txt

OPTIONS   =  -M0

all: $(OUT_FILES)

uterus.txt: uterus.rxt
	$(IAD_EXECUTABLE) -g 0.9 uterus.rxt
	
sevick-b.txt: sevick-b.rxt
	$(IAD_EXECUTABLE) -g 0.9 sevick-b.rxt

sevick-A.txt: sevick-A.rxt
	$(IAD_EXECUTABLE) -g 0.9 sevick-A.rxt

ville1.txt: ville1.rxt
	$(IAD_EXECUTABLE) -a 0 ville1.rxt

basic-E.txt: basic-E.rxt
	$(IAD_EXECUTABLE) -g 0.9 basic-E.rxt

$(OUT_FILES): %.txt: %.rxt
	$(IAD_EXECUTABLE) $(OPTIONS) $<

clean:
	rm -f *.txt

.PHONY : clean
