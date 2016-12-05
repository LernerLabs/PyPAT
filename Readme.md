#Installation

New:

Install locally with

```bash
pip install -e .
```

Old:

The executable scripts are located in the "drivers" subdirectory, and
the libraries are located in the "pypat" directory. Inside each of the
scripts, you will need to edit the lines containing
"PYPAT_CODE_DIRECTORY" so that they point to the directory in which
the PyPAT package is installed.

#Dependencies

Python 2.5 or greater is required.
Your Python installation must have the following installed:

 - matplotlib 0.9.0
 - numpy 1.0.1
 - PyMOL 1.0

And you must have following installed on your system:

 - gnuplot
 - ImageMagick 6.3.2
 - PyMOL 1.0

#Running the Tools

The executable scripts are found in the "drivers" subdirectory. You
may wish to place this directory in your PATH. All scripts support the
--help command line option. Documentation for individual tools
follows.

TODO: get this from the wiki
