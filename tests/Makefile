#
# CLI fixture batch tests. Command-line option coverage lives in tests/cli.
#

IAD_EXECUTABLE ?= ../iad

all:
	IAD_EXECUTABLE=$(IAD_EXECUTABLE) cli/run_cli_tests.sh batch

clean:
	rm -f *.txt *.svg *.ipynb *.grid
	rm -rf *notebook

.PHONY : all clean
