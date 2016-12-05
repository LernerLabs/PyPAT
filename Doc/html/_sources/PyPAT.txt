PyPAT
=====

PyPAT (Python-based Protein Analysis Tools) is a collection of tools
that build upon the ptraj module of AMBER and the PyMOL visualization
package to aid in the analysis of protein structures and molecular
dynamics trajectories. They allow for the evaluation of the
convergence of trajectories, as well as the examination of correlated
dynamics, hydrogen bonds, and bridging-water molecules throughout a
trajectory. Our tools are written in Python and released under an
open-source license.

This document is intended to explain the usage of the PyPAT tools. For
an understanding of the theory behind, and implementation of, these
tools, users are refered to the paper:

Lerner, M.G., Spronk, S.A., and Carlson, H.A. (2008) PyPAT: a
Python-based toolset to aid in the analysis of protein structures and
trajectories. *submitted to Bioinformatics*.


System Requirements
===================

Our tools are primarily meant for use with the AMBER suite of
programs, and they have been extensively tested on both OS X and
Linux. Many may be installed on Windows systems, but this is generally
unsupported.

Required software includes

 - gnuplot_ 4.2 (www.gnuplot.info)

 - ImageMagick_ 6.4.3 (www.imagemagick.org),

 - matplotlib_ 0.98.3 (matplot-lib.sourceforge.net)

 - numpy_ 1.1.1 (numpy.scipy.org)

 - PyMOL_ 1.1 (www.pymol.org)

 - Python_ 2.6 (www.python.org)

all of which are freely available and open-source.

Installation
============

Installation on Linux and OS X follows standard Python proctocols::

  prompt$ python setup.py build
  prompt$ python setup.py install

Users may easily install into a local directory with the ``--prefix`` option::

  prompt$ python setup.py install --prefix=/my/home/dir/software

Documentation
=============

Documentation is provided in `ReStructured Text`_ formatted text files
in the top-level directory. Sphinx_ has been used to generate html and
pdf documentation from these text files. The html and pdf files live
in the ``Doc`` subdirectory. If the need arises, documentation may be
rebuilt::

  prompt$ make clean
  prompt$ make html
  prompt$ make latex
  prompt$ cd .build/latex
  prompt$ make all-pdf
  prompt$ cd ../..
  prompt$ rm -rf Doc/html
  prompt$ mv .build/html Doc
  prompt$ mv .build/latex/PyPAT.pdf Doc

.. _`ReStructured Text`: http://docutils.sourceforge.net/docs/user/rst/quickref.html
.. _Sphinx: http://sphinx.pocoo.org/
.. _Python: http://www.python.org
.. _matplotlib: http://matplot-lib.sourceforge.net
.. _numpy: http://numpy.scipy.org
.. _PyMOL: http://www.pymol.org
.. _gnuplot: http://www.gnuplot.info
.. _ImageMagick: http://www.imagemagick.org
