OVERVIEW
========

`ad` and `iad` are command-line utilities for radiative transport in layered geometries.  `ad` is short for van de Hulst's adding-doubling method for finding the total reflection and transmission through a layered medium.  `iad` is short for inverse adding-doubling that solves the inverse problem of finding the intrinsic optical properties from the total reflection and transmission.

INSTALLATION
============

On Linux and macOS, you should be able to just type

    make

to create excutable versions of the `ad` and `iad` programs.  

> Despite being written in standard C, cross-platform compiling can be a nuisance because I wrote the program using the literate [cweb](http://literateprogramming.com/cweb_download.html) techniques of Knuth.  This means the original source files `*.w` need to be converted to to `*.c` files using a program called `ctangle`. Since I don't want everyone else to be forced to wear the `cweb` hairshirt, I have included the `.c` and `.h` files, but sometimes I screw up.  The ctangle code is available from [Stanford](https://www-cs-faculty.stanford.edu/~knuth/cweb.html).


To verify that the program is working properly, try

    make test

to initiate a series of command-line tests.  All the tests should pass (there'll
be a '*' at the end of each line).  Alternatively, you can try these from the
command-line

    ./ad -a 0.8 -b 1.0 -g 0.5

    ./iad -r 0.09636 -t 0.67004 -d 0.5

    ./iad -M 0 -q 4 test/vio-A

The last of these processes the file `test/vio-A.rxt` to create `test/vio-A.txt`.  To install the program

    sudo make install

and the binary executable files (`iad` and `ad`) will be created and placed in 
`$(DESTDIR)/bin`.
