#
#  Makefile by Scott Prahl in Nov 1998
#

IAD_EXECUTABLE ?= ../iad

TST_FILES = basic_A.rxt     basic_B.rxt    basic_C.rxt    basic_D.rxt    basic_E.rxt   \
            double.rxt      example2.rxt   fairway_A.rxt  fairway_B.rxt  fairway_C.rxt \
            fairway_D.rxt   fairway_E.rxt  il_A.rxt       il_B.rxt       il_C.rxt      \
            ink_A.rxt       ink_B.rxt      ink_C.rxt      kenlee_A.rxt   kenlee_B.rxt  \
            kenlee_C.rxt    newton.rxt     royston1.rxt   royston2.rxt   royston3_A.rxt\
            royston3_B.rxt  royston3_C.rxt royston3_D.rxt royston3_E.rxt royston9_A.rxt\
            royston9_B.rxt  royston9_C.rxt royston9_D.rxt sample_A.rxt   sample_B.rxt  \
            sample_C.rxt    sample_D.rxt   sample_E.rxt   sample_F.rxt   sample_G.rxt  \
            sevick_A.rxt    sevick_B.rxt   terse_A.rxt    terse_B.rxt    tio2_vis.rxt  \
            uterus.rxt      ville1.rxt     vio_A.rxt      vio_B.rxt      x_bad_data.rxt

OUT_FILES = basic_A.txt    basic_B.txt    basic_C.txt    basic_D.txt    basic_E.txt    \
            double.txt     example2.txt   fairway_B.txt  fairway_C.txt  fairway_D.txt  \
            il_A.txt       ink_A.txt      ink_B.txt      ink_C.txt      kenlee_A.txt   \
            kenlee_B.txt   kenlee_C.txt   newton.txt     royston1.txt   royston2.txt   \
            royston3_A.txt royston3_B.txt royston3_C.txt royston3_D.txt royston3_E.txt \
            royston9_A.txt royston9_B.txt royston9_C.txt  \
            sample_A.txt   sample_B.txt   sample_C.txt   sample_D.txt   sample_E.txt   \
            sample_F.txt   sample_G.txt   terse_A.txt    terse_B.txt    vio_A.txt      \
            vio_B.txt      x_bad_data.txt il_B.txt       il_C.txt

SPECIAL_FILES = ville1.txt 

NO_MC_FILES = royston1.txt tio2_vis.txt royston9_D.txt

NO_RC_FILES = fairway_A.txt fairway_E.txt

OPTIONS   =  

all: $(OUT_FILES) $(SPECIAL_FILES) $(NO_MC_FILES) $(NO_RC_FILES)

$(OUT_FILES): %.txt: %.rxt
	$(IAD_EXECUTABLE) $(OPTIONS) $<

ville1.txt: ville1.rxt
	$(IAD_EXECUTABLE) -a 0 $<

# These files all have too little reflected light collected
fairway_A.txt: fairway_A.rxt
	$(IAD_EXECUTABLE) -c 0 -M 0 $<

fairway_E.txt: fairway_E.rxt
	$(IAD_EXECUTABLE) -c 0 -M 0 $<

# These files fail when lost light is calculated
tio2_vis.txt: tio2_vis.rxt
	$(IAD_EXECUTABLE) -M 0 $<

royston9_D.txt: royston9_D.rxt
	$(IAD_EXECUTABLE) -e 0.005 $<

clean:
	rm -f *.txt *.svg *.ipynb *.grid
	rm -rf *notebook

.PHONY : clean
