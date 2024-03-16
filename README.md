# Inverse Adding-Doubling

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/scottprahl/iad?label=latest)
[![MIT License](https://img.shields.io/badge/MIT-license-yellow.svg)](https://github.com/scottprahl/miepython/blob/master/LICENSE.txt)
[![GitHub Actions](https://github.com/scottprahl/iad/actions/workflows/make.yml/badge.svg)](https://github.com/scottprahl/iad/actions/workflows/make.yml)
[![DOI](https://zenodo.org/badge/102147394.svg)](https://zenodo.org/badge/latestdoi/102147394)

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

### Jupyter support

As of March 2024, there is now a python command-line script `iadplus` that will analyze an `.rxt` input file and graph the results.  Everything is assembled into a Jupyter notebook for convenience.  You may need to install some python modules to be able to use `iadplus`

    iadplus -options '-i 8 -X ' file.rxt

will produce `file.txt` as well as bunch of `.svg` files for the Jupyter notebook
`file.ipynb`

## Author

Scott Prahl

http://omlc.org/~prahl
