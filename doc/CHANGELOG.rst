Changelog
=========

v3.12.0 (unreleased)
-------------------
*   add continuous building (tvercat)
*   improve cweave/ctwill processing (ascherer)
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
