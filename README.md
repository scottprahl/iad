# Inverse Adding-Doubling

** VERSION 3.10.0 **

** Scott Prahl **

** Oregon Tech **

## OVERVIEW

Inverse Adding-Doubling (IAD) determines the intrinsic optical properties of a scattering and absoption sample using measurements of the total reflection and transmission.  Basically, optical properties are repeatedly guessed until the calculated reflection and transmission match the measured values.

This package provides two executables `ad` and `iad`.  The first does a forward adding-doubling calculation (i.e., given the albedo, optical thickness, and anisotropy it returns the total reflection and transmission).  The second
does the reverse.

This program (Prahl et al., *Applied Optics*, **32**, 559-568, 1993) uses the Adding-Doubling method of
van de Hulst (*Multiple Light Scattering*, Academic Press, 1978).  The Adding-Doubling method was extended to account for index of refraction changes at boundaries as well as corrections that must accompany integrating sphere experiments.

Finally, integrating spheres do not always collect all the light that exits from the front or back surface of a sample.  Since this is impossible to account for the 1D adding-doubling technique, a Monte Carlo simulation is included in the inverse calculation.


## INSTALLATION

In principle you should be able to just type

    make

to create excutable versions of the `ad` and `iad` programs.  To quickly test
if these are working try

    ./ad -a 0.8 -b 1.0 -g 0.5

    ./iad -r 0.09636 -t 0.67004 -d 0.5

    ./iad -M 0 -q 4 test/vio-A

The last command will process the file `vio-A.rxt` to produce `vio-A.abg`.  To install

    make install

and the binary executable files (`iad` and `ad`) will be created and placed in a default location of
`/usr/local/bin/`.


### Shared library support.  

Edit the Makefile to select the right type of shared library for your platform

	make install-lib

### Mathematica support.  

If you have Mathematica, then (and only if you have installed the right
tools and edited the Makefile for your platform) and have the libraries installed then you should be able to type

	make mma
	make install mma

and then load the iad module and then type 

    Plot[UR1[0.5,1.0,g], {g,0.5,0.8}]

in Mathematica to get a graph.  Very cool.

### Python support

Once you have installed the shared library (`.dylib` under MacOS) or (`.so` under linux)
then you can install python bindings

    pip install iadpython

then in Jupyter 

    import numpy as np
    import matplotlib.pyplot as plt
    import iadpython as iad
    
    g = np.linspace(0.5,0.8,50)
    plt.plot(g, iad.rt(1,1,0.5,1.0,g))
    plt.show()
    

## Author

Scott Prahl

`scott.prahl@oit.edu`

http://omlc.org/~prahl
