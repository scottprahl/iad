# Inverse Adding-Doubling

[![GitHub Code](https://img.shields.io/badge/github-code-green.svg)](https://github.com/scottprahl/iad)
[![MIT License](https://img.shields.io/badge/MIT-license-yellow.svg)](https://github.com/scottprahl/miepython/blob/master/LICENSE.txt)
[![GitHub Actions](https://github.com/scottprahl/iad/actions/workflows/make.yml/badge.svg)](https://github.com/scottprahl/iad/actions/workflows/make.yml)
[![version](https://img.shields.io/github/v/release/scottprahl/iad.git)](https://img.shields.io/github/v/tag/scottprahl/iad)


by Scott Prahl

May 2023

## OVERVIEW

Inverse Adding-Doubling is a command-line program that determines the intrinsic optical properties of a flat scattering and absoption sample using measurements of the total reflection and transmission.  Basically, optical properties are repeatedly guessed until the calculated reflection and transmission match the measured values.

This package provides two executables `ad` and `iad`.  The first does a forward adding-doubling calculation (i.e., given the albedo, optical thickness, and anisotropy it returns the total reflection and transmission).  The second
does the reverse.

This program [Prahl et al., *Applied Optics*, **32**, 559-568, 1993](https://omlc.org/~prahl/pubs/pdfx/prahl93a.pdf) uses the Adding-Doubling method of [van de Hulst *Multiple Light Scattering*, Academic Press, 1978](https://www.amazon.com/Multiple-Light-Scattering-Formulas-Applications-ebook/dp/B01D4CMF80).  I extended the Adding-Doubling method to account for Fresnel reflection at boundaries as well as corrections that must accompany integrating sphere experiments.

Finally, integrating spheres do not always collect all the light that exits from the front or back surface of a sample.  Since this is impossible to account for the 1D adding-doubling technique, a Monte Carlo simulation is included in the inverse calculation.

Details about using the program are documented in the accompanying [manual](/doc/manual.pdf).

## INSTALLATION

In principle, in a unix environment you should be able to just type

    make install

to create and install executable versions of the `ad` and `iad` programs.  See
[INSTALL.md](/INSTALL.md) for more details. Then run

    iad test/basic-A.rxt

to translate the reflection and transmission measurements to optical properties in the generated file `test/basic-A.txt`

> For Windows, there are executable binaries `ad.exe` and `iad.exe` compiled using [MinGW-w64](https://mingw-w64.org/doku.php).  These apps can be run using the `Command Prompt` application `cmd.exe`.  These binaries are packaged in a separate `iad-win` distributions on [github](https://github.com/scottprahl/iad/releases) or [omlc](https://omlc.org/software/iad/).

### Changelog

This is located in the `doc/` directory.

### Python support

Once you have installed the shared library (`.dylib` under macOS) or (`.so` under linux) then you can install python bindings

    pip install iadpython

then in Jupyter 

    import numpy as np
    import matplotlib.pyplot as plt
    import iadpython as iad
    
    g = np.linspace(0.5,0.8,50)
    plt.plot(g, iad.rt(1,1,0.5,1.0,g))
    plt.show()

### Shared library support.  

Edit the `Makefile` to select the right type of shared library for your platform

    make install-lib

### Mathematica support.  

If you have Mathematica, then (and only if you have installed the right
tools and edited the Makefile for your platform) and have the libraries installed then you should be able to type

    make mma
	make install mma

and then load the iad module and then type 

    Plot[UR1[0.5,1.0,g], {g,0.5,0.8}]

in Mathematica to get a graph.  Very cool.

## Author

Scott Prahl

http://omlc.org/~prahl
