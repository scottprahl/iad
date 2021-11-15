#
#  Makefile by Scott Prahl in Nov 2010
#
#  Requires the IAD library files to be installed 
#  Specifically a header and library are needed
#            libiad.a     and libiad.h
#

VERSION      = 8.0
SYS          = MacOSX

CC           = /usr/local/Cellar/gcc
INSTALL_DIR  = /Users/prahl/Library/Mathematica/Applications/Optics
CADDSDIR     = /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/CompilerAdditions

INC_DIR      = ${CADDSDIR}
LIB_DIR      = ${CADDSDIR}
EXTRA_LIBS   = -lstdc++ -framework Foundation -liad 
MPREP        = ${CADDSDIR}/mprep

all : AD

.c.o :
	${CC} -c -I${INC_DIR} -I/usr/local/include $<

ADtm.c : AD.tm
	${MPREP} $? -o $@

AD : ADtm.o
	${CC} -I${INC_DIR} ADtm.o -L${LIB_DIR} -lMLi4 ${EXTRA_LIBS} -o $@

install:
	mkdir -p ${INSTALL_DIR}
	cp AD.m ${INSTALL_DIR}
	mkdir -p ${INSTALL_DIR}/External
	cp AD ${INSTALL_DIR}/External

clean:
	rm -f ADtm.o AD.o ADtm.c AD
	
.PHONY : clean install lib
