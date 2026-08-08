Changelog
=========

v4.0.0 (29 Jul 2026)
---------------------

Backwards-incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*   remove the faulty power-law forms of ``-F`` (``-F 'P 500 1.0 -1.3'`` and
    ``-F 'R 500 1.0 -1.3'``); the other uses of ``-F`` are unchanged
*   test fixtures moved from ``test/`` to ``tests/rxt/{0,1,2}_sphere``

Inversion
~~~~~~~~~
*   replace the precomputed grid with an adaptive grid (``iad_agrid``) for initial guesses
*   support ``FIND_BaG`` and ``FIND_BsG`` searches with the adaptive grid
*   add new L2_SCALED deviation metric and make it the default
*   add bounded Nelder-Mead amoeba with explicit lower/upper bounds
*   initialize the simplex in physical space using scipy-style 5% steps
*   add Monte Carlo simplex hot-start so re-inversions reuse previous solution
*   reduce MC update factor from 0.8 to 0.3
*   no longer try to rescue a failed convergence with Monte Carlo
*   accept a valid boundary-clamped solution even after the simplex hits its iteration cap
*   improved single sphere MT calc when unscattered T exits
*   reset defaults for each data point
*   use a non-deterministic Monte Carlo seed each run

Command line and output
~~~~~~~~~~~~~~~~~~~~~~~
*   ``iad *`` processes every ``.rxt`` file matched by the glob and ignores the rest
*   accept several input files in one invocation and print each filename as it is processed
*   ``iad -z`` now prints intrinsic, derived, sphere, and calculated R/T quantities
*   allow incidence angle to appear in a column of ``.rxt`` files
*   when ``-w`` is given, use it for ``-W`` as well
*   tag output header with ``iad`` so the writing program is identifiable
*   clearer messages when an inversion fails
*   catch low M_R measurements early
*   suppress spurious "max R" error with double spheres

Fixes
~~~~~
*   fix ``-o outfile`` placed after the input filename (POSIX getopt on macOS/BSD)
*   fix missing break in ``-x`` option parsing
*   flag ``M_R + M_T > 1`` with ``!`` when no spheres are used, instead of
    silently reporting a bogus fit (e.g. ``iad -r 0.4 -t 0.7``)
*   explain failures for measurements given on the command line; nothing was
    printed before.  A single ``-r``/``-t``/``-u`` measurement now names its
    one error, while input files still get the full legend
*   give every error code its own message instead of silently printing
    nothing for the sphere, phase-function, and file errors
*   show ``AD`` rather than ``IAD`` above the iterations for ``-x 8``
*   reallocate redistribution-function cache when quadrature-point count changes
*   fix ``mc_lost`` build
*   single-sphere command-line fixes
*   improve M_T calculation
*   improve lost-light and sphere debug output

Documentation
~~~~~~~~~~~~~
*   rewrite the README and convert it from Markdown to ``README.rst``
*   document every command-line option in the manual
*   describe the current inversion algorithm in the manual
*   ``cweave`` now produces warning-free docs; fix ``iad_src.pdf`` generation

Build, source, and tests
~~~~~~~~~~~~~~~~~~~~~~~~
*   fix all compiler warnings
*   convert ``mc_lost.c`` to a CWEB ``mc_lost.w`` source
*   rewrite C-style comments in the ``.w`` files in CWEB style
*   rename source files so that no filename contains a hyphen
*   move command-line tests out of the Makefile into shell scripts under ``tests/cli``
*   generate ``src/version.c`` with make and retire ``version.pl``
*   ``make dist`` and ``make windist`` follow the new test layout and ship
    ``iadplus``, ``iadsum``, and ``CITATION.cff``
*   remove dead code
*   tabs-to-spaces and formatting passes across CWEB sources

Python helpers
~~~~~~~~~~~~~~
*   add ``--html`` output option and docstrings to ``iadplus``
*   improve the notebooks and plots that ``iadplus`` creates
*   report a clear error when the ``iad`` executable cannot be found
*   improve ``iadsum``

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
---------------------
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
---------------------
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
---------------------
*   produce 64-bit windows executable since
*   32-bit triggered false positive virus detection

v3.14.1 (30 Jan 2024)
---------------------
*   no longer toss correct solution in some cases
*   only calculate redistribution matrix when needed
*   improve debug comments
*   start stripping tabs from cweb files

v3.14.0 (25 Jan 2024)
---------------------
*   fix handling of slides (@anishabahl)
*   fix github build
*   improve Makefile
*   warn on bad sphere wall reflectivity

v3.13.2 (24 Jan 2024)
---------------------
*   fix port size normalization (@jgroehl)
*   update copyright year

v3.13.1 (24 Jan 2024)
---------------------
*   left debugging statements in

v3.13.0 (24 Jan 2024)
---------------------
*   add -1 feature for parameters in .rxt files

v3.12.1 (26 May 2023)
---------------------
*   bump version to get zenodo links correct

v3.12.0 (26 May 2023)
----------------------
*   add continuous building (@tvercat)
*   improve cweave/ctwill processing (@ascherer)
*   add CITATION.cff to base level of repository
*   add DOI for citation purposes
*   added badges to README page (whee!)

v3.11.6 (18 Nov 2021)
---------------------
*   fix initialization for couple of corner cases (finding just g)
*   found while adding tests to iadpython.

v3.11.5 (7 Nov 2020)
--------------------
*   fix initialization problem when using ad_layers

v3.11.4 (15 Oct 2019)
---------------------
*   solve compilation problem on Raspberry Pi by adding -fsigned-char complier option
*   touch .c and .h files in Makefile to avoid needing ctangle

v3.11.3 (21 Aug 2019)
---------------------
*   improve an error message when using -F
*   add command-line option to specify search explicitly
*   improve help message

v3.11.2 (29 Mar 2019)
---------------------
*   Add separate License file
*   Make copyright notices consistent
*   Add some basic hints to the README.md for Windows users
*   Update the doc/CHANGELOG

v3.11.1 (28 Mar 2019)
---------------------
*   The main change in this release is that windows executables can now be built with MinGW-w64 and tested under Wine.

v3.10.3 (3 July 2018)
---------------------
*   This release mostly improves packaging so that everything compiles cleanly on MacOS X and linux.
*   Improved tests and fixed a few minor bugs in the frameworks
*   Improved information presented during debugging.

v3.10.2 (2 Nov 2017)
--------------------
*   This version adds header files needed to install libiad that formerly needed to be generated with ctangle.

v3.10.1 (2 Nov 2017)
--------------------
*   This version now includes .c and .h files that are generated using the ctangle program. The program should build cleanly on unix/macos platforms.
