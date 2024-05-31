# Inverse Adding-Doubling

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/scottprahl/iad?label=latest)
[![MIT License](https://img.shields.io/badge/MIT-license-yellow.svg)](https://github.com/scottprahl/miepython/blob/main/LICENSE.txt)
[![GitHub Actions](https://github.com/scottprahl/iad/actions/workflows/make.yml/badge.svg)](https://github.com/scottprahl/iad/actions/workflows/make.yml)
[![DOI](https://zenodo.org/badge/102147394.svg)](https://zenodo.org/badge/latestdoi/102147394)

by Scott Prahl

April 2024

## Overview

Inverse Adding-Doubling is a command-line program that determines the intrinsic optical properties of a flat scattering and absoption sample using measurements of the total reflection and transmission.  Basically, optical properties are repeatedly guessed until the calculated reflection and transmission match the measured values. The program accounts
for realistic measurement constraints associated with integrating sphere effects.

This package provides two executables `ad` and `iad`.  The first does a forward adding-doubling calculation (i.e., given the albedo, optical thickness, and anisotropy it returns the total reflection and transmission).  The second does the reverse.

This program [Prahl et al., *Applied Optics*, **32**, 559-568, 1993](https://omlc.org/~prahl/pubs/pdfx/prahl93a.pdf) uses the Adding-Doubling method of [van de Hulst *Multiple Light Scattering*, Academic Press, 1978](https://www.amazon.com/Multiple-Light-Scattering-Formulas-Applications-ebook/dp/B01D4CMF80).  I extended the Adding-Doubling method to account for Fresnel reflection at boundaries as well as corrections that must accompany integrating sphere experiments.

Finally, integrating spheres do not always collect all the light that exits from the front or back surface of a sample.  Since this is impossible to account for the 1D adding-doubling technique, a Monte Carlo simulation is embedded in the inverse calculation.

Details about using the program are documented in the accompanying [manual](/doc/manual.pdf).

## Installation

Obtaining executables from source entails unzipping and running `make` (See
[INSTALL.md](/INSTALL.md) for more details.) Windows executables are available
in the `iad-win-3-16-2.zip` releases.

## Usage

### Inverting a single measurement on the command line

To find the optical properties of a 1mm thick sample with 40% total reflectance and 10% total transmission

```bash
    ./iad -r 0.4 -t 0.1 -d 1
```

which will generate output that ends with

```
#       Measured   M_R     Measured   M_T    Estimated Estimated Estimated
##wave     M_R     fit        M_T     fit       mu_a     mu_s'       g    
# [nm]    [---]    [---]     [---]    [---]     1/mm     1/mm      [---]  
     1   0.4000   0.4000    0.1000   0.0999    0.4673   3.9317    0.0000    *
```

The process can be reversed by doing a forward calculation with the `-z` option.  Here `-A 0.4673` sets the absorption coefficient, and `-j 3.9317` sets the reduced scattering coefficient).  Specifying the thickness `-d 1` is necessary because the default
thickness is infinite.


```bash
    ./iad -z -A 0.4673 -j 3.9317 -d 1
```

produces output that ends with

```
#       Measured    M_R    Measured    M_T    Estimated Estimated Estimated
##wave     M_R      fit       M_T      fit      mu_a      mu_s'       g    
# [nm]    [---]    [---]     [---]    [---]     1/mm      1/mm      [---]  
     0   0.0000   0.4000    0.0000   0.0999    0.4673    3.9317     0.0000    * 
```

> For Windows, there are executable binaries `ad.exe` and `iad.exe` compiled using [MinGW-w64](https://mingw-w64.org/doku.php).  These apps can be run using the `Command Prompt` application `cmd.exe`.  These binaries are packaged in the separate `iad-win` distributions on [github](https://github.com/scottprahl/iad/releases) or [omlc](https://omlc.org/software/iad/).

### Inverting spectral data

Often one wants the optical properties over an entire spectrum.  A good example was
recently provide by @anishabahl.  This measurement was made with a spectrophotometer 
equipped with a dual beam integrating sphere.  The [input data](https://github.com/scottprahl/iad/blob/main/docs/phantom-with-no-slides.rxt) looks includes the total reflection and transmission as well a header that describes the
experiment.  The reflection and transmission data look like this

![r and t graph](https://github.com/scottprahl/iad/blob/main/docs/phantom-with-no-slides-RTU.svg)

The input file is processed with

```bash
    iad -X -i 8 -g 0.9 phantom-with-no-slides.rxt
```

The option `i 8` indicates that light is incident on the sample at an angle of 8°, `-X` indicates that a sphere with dual beams was used, and `-g 0.9` indicates  the default
scattering anisotropy.  Note that in the PDMS file, the refractive index of the
sample changes with each data point

This command will produce an [output file](https://github.com/scottprahl/iad/blob/main/docs/phantom-with-no-slides.txt) that when plotted looks like

![calculated mua](https://github.com/scottprahl/iad/blob/main/docs/phantom-with-no-slides-mua.svg)

which matches the 
intrinsic absorption of Wacker PDMS found in Cai's 2008 [Disseration](https://dx.doi.org/10.17877/DE290R-8242)

![calculated mus](https://github.com/scottprahl/iad/blob/main/docs/phantom-with-no-slides-mus.svg)

This is an excellent set of measurements because there is
almost no influence of the absorption coefficient on the scattering coefficient.  
 
### Jupyter support

As of March 2024, there is now a python command-line script `iadplus` that will analyze an `.rxt` input file and graph the results.  Everything is assembled into a Jupyter notebook for convenience.  You may need to install some python modules to be able to use `iadplus`

    iadplus -options '-i 8 -X ' file.rxt

will run `iad` and produce `file-notebook/file.ipynb` along with a bunch of `.svg` files used in the Jupyter notebook  `file.ipynb`

## Author

Scott Prahl

http://omlc.org/~prahl
