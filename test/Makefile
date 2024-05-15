#
#  Makefile by Scott Prahl in Nov 1998
#

IAD_EXECUTABLE ?= ../iad

TST_FILES = basic-A.rxt     basic-B.rxt    basic-C.rxt    basic-D.rxt    basic-E.rxt   \
            double.rxt      example2.rxt   fairway-A.rxt  fairway-B.rxt  fairway-C.rxt \
            fairway-D.rxt   fairway-E.rxt  il-A.rxt       il-B.rxt       il-C.rxt      \
            ink-A.rxt       ink-B.rxt      ink-C.rxt      kenlee-A.rxt   kenlee-B.rxt  \
            kenlee-C.rxt    newton.rxt     royston1.rxt   royston2.rxt   royston3-A.rxt\
            royston3-B.rxt  royston3-C.rxt royston3-D.rxt royston3-E.rxt royston9-A.rxt\
            royston9-B.rxt  royston9-C.rxt royston9-D.rxt sample-A.rxt   sample-B.rxt  \
            sample-C.rxt    sample-D.rxt   sample-E.rxt   sample-F.rxt   sample-G.rxt  \
            sevick-A.rxt    sevick-b.rxt   terse-A.rxt    terse-B.rxt    tio2_vis.rxt  \
            uterus.rxt      ville1.rxt     vio-A.rxt      vio-B.rxt      x_bad_data.rxt

OUT_FILES = basic-A.txt    basic-B.txt    basic-C.txt    basic-D.txt    basic-E.txt    \
            double.txt     example2.txt   fairway-B.txt  fairway-C.txt  fairway-D.txt  \
            il-A.txt       ink-A.txt      ink-B.txt      ink-C.txt      kenlee-A.txt   \
            kenlee-B.txt   kenlee-C.txt   newton.txt     royston1.txt   royston2.txt   \
            royston3-A.txt royston3-B.txt royston3-C.txt royston3-D.txt royston3-E.txt \
            royston9-A.txt royston9-B.txt royston9-C.txt  \
            sample-A.txt   sample-B.txt   sample-C.txt   sample-D.txt   sample-E.txt   \
            sample-F.txt   sample-G.txt   terse-A.txt    terse-B.txt    vio-A.txt      \
            vio-B.txt      x_bad_data.txt il-B.txt       il-C.txt

SPECIAL_FILES = ville1.txt 

NO_MC_FILES = royston1.txt tio2_vis.txt royston9-D.txt

NO_RC_FILES = fairway-A.txt fairway-E.txt

OPTIONS   =  

all: $(OUT_FILES) $(SPECIAL_FILES) $(NO_MC_FILES) $(NO_RC_FILES)

$(OUT_FILES): %.txt: %.rxt
	$(IAD_EXECUTABLE) $(OPTIONS) $<

ville1.txt: ville1.rxt
	$(IAD_EXECUTABLE) -a 0 $<

# These files all have too little reflected light collected
fairway-A.txt: fairway-A.rxt
	$(IAD_EXECUTABLE) -c 0 -M 0 $<

fairway-E.txt: fairway-E.rxt
	$(IAD_EXECUTABLE) -c 0 -M 0 $<

# These files fail when lost light is calculated
tio2_vis.txt: tio2_vis.rxt
	$(IAD_EXECUTABLE) -M 0 $<

royston9-D.txt: royston9-D.rxt
	$(IAD_EXECUTABLE) -e 0.005 $<

clean:
	rm -f *.txt *.svg *.ipynb *.grid
	rm -rf *notebook

.PHONY : clean
