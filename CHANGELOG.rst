Changelog
=========

v3.16.4 (unreleased)
---------------------
*   improved single sphere MT calc when unscattered T leaves

v3.16.3 (15 May 2024)
---------------------
*   fixed bug that discarded info from "g" input column
*   improved readme

v3.16.2 (23 Apr 2024)
---------------------
*   fixed bug in mc_lost affecting lost light estimation
*   force one MC run when port sizes are present
*   revised manual
*   add wavelength constraints -l '500 600' to limit processing
*   allow reduced scattering to be specified with -j
*   better checks for out of date redistribution matrix
*   exit properly with bad header
*   reduce max iterations to 100
*   include mc_test.c in distribution
*   show line number for bad .rxt entry
*   invalidate h calc when angles change
*   build targets for test programs
*   better debug -x 1 output
*   improved iadsum
*   improved iadplus 
*   made mc_lost as standalong executable
*   made forward calculation work with no spheres
*   improved g spacing in grid
*   better support for m.u measurements

v3.16.1 (25 Mar 2024)
--------------------
*   clarify and revise single sphere effects
*   avoid MC for failed 1 parameter searches
*   change 'empty' -> 'entrance' or 'empty' -> 'third' as appropriate
*   enable constraints in rxt files
*   fix command line regression of constraints
*   improve -x 2 debugging for grid generation
*   improve Valid_Grid
*   add -J to generate grid for plotting
*   improvements to iadplus

v3.16.0 (16 Mar 2024)
--------------------
*   flexible columns input files!
*   add -w and -W command line options
*   include MC lost light in -z calculations
*   include MC lost light in 1 parameter searches
*   better support for overwriting r&t file data from command-line
*   first version of python wrapper iadplus
*   add more simple checks for bad r & t values

v3.15.1 (4 Mar 2024)
--------------------
*   add -x 1 support back for debugging
*   add -L 633 to set wavelength to 633nm
*   fix building docs/iad_src.pdf
*   add separate mc_test and mc_lost_test
*   improve -x debugging
*   allow more constraint situations
*   remove unused code

v3.15.0 (2 Mar 2024)
--------------------
*   support for spheres with no baffles
*   better help message
*   works better with extreme anisotropies
*   clarify MT calculations
*   improve header in output file

v3.14.5 (5 Feb 2024)
--------------------
*   improve command line parameter handling
*   fix URU calc for oblique angles
*   use exit(EXIT_SUCCESSS) and exit(EXIT_FAILURE) consistently
*   detab more files
*   add mctest target

v3.14.4 (1 Feb 2024)
--------------------
*   fix lost diffuse light code
*   dramatically improve the look of generated .c and .h files
*   remove Mathematica support
*   Makefile cleanup
*   add test for lost light calculation
*   use POSIX getopt
*   use M_PI instead of number
*   add test code for lost light estimation
*   fix oblique test code

v3.14.3 (31 Jan 2024)
--------------------
*   produce 64-bit windows executable since
*   32-bit triggered false positive virus detection

v3.14.1 (30 Jan 2024)
--------------------
*   no longer toss correct solution in some cases
*   only calculate redistribution matrix when needed
*   improve debug comments
*   start stripping tabs from cweb files

v3.14.0 (25 Jan 2024)
--------------------
*   fix handling of slides (@anishabahl)
*   fix github build
*   improve Makefile
*   warn on bad sphere wall reflectivity

v3.13.2 (24 Jan 2024)
--------------------
*   fix port size normalization (@jgroehl)
*   update copyright year

v3.13.1 (24 Jan 2024)
--------------------
*   left debugging statements in

v3.13.0 (24 Jan 2024)
--------------------
*   add -1 feature for parameters in .rxt files

v3.12.1 (26 May 2023)
---------------------
*   bump version to get zenodo links correct

v3.12.0
-------------------
*   add continuous building (@tvercat)
*   improve cweave/ctwill processing (@ascherer)
*   add CITATION.cff to base level of repository
*   add DOI for citation purposes
*   added badges to README page (whee!)

v3.11.6
-------------------
*   fix initialization for couple of corner cases (finding just g)
*   found while adding tests to iadpython.

v3.11.5
-------------------
*   fix initialization problem when using ad_layers

v3.11.4
-------------------
*   solve compilation problem on Raspberry Pi by adding -fsigned-char complier option
*   touch .c and .h files in Makefile to avoid needing ctangle

v3.11.3
-------------------
*   improve an error message when using -F
*   add command-line option to specify search explicitly
*   improve help message

v3.11.2
-------------------
*   Add separate License file
*   Make copyright notices consistent
*   Add some basic hints to the README.md for Windows users
*   Update the doc/CHANGELOG

v3.11.1
-------------------
*   The main change in this release is that windows executables can now be built with MinGW-w64 and tested under Wine.

v10.3.3
-------------------
*   This release mostly improves packaging so that everything compiles cleanly on MacOS X and linux.
*   Improved tests and fixed a few minor bugs in the frameworks
*   Improved information presented during debugging.

v10.3.2
-------------------
*   This version adds header files needed to install libiad that formerly needed to be generated with ctangle.

v3.10.1
-------------------
*   This version now includes .c and .h files that are generated using the ctangle program. The program should build cleanly on unix/macos platforms.
