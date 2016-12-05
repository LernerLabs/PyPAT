#!/usr/bin/env python

from distutils.core import setup

setup(name='PyPAT',
      version='1.0',
      description='The PyPAT package',
      author='Michael G. Lerner, Steven A. Spronk, Heather A. Carlson',
      author_email='mglerner@gmail.com',
      url='http://github.com/mglerner/PyPAT',
##      packages=['pypat', 'pypat.hbond', 'pypat.runningptraj'],
      # The magic thing about scripts is that "if the first line of
      # the script starts with #! and contains the word ``python'',
      # the Distutils will adjust the first line to refer to the
      # current interpreter location."
      # (http://docs.python.org/dist/node11.html)
##      scripts=['drivers/collect_water_bridges.py',
##               'drivers/convert_to_numpy_format.py',
##               'drivers/display_bridging_interactions.py',
##               'drivers/do_correlated_md_analysis.py',
##               'drivers/make_correlated_dynamics_plots.py',
##               'drivers/make_movies.py',
##               'drivers/parse_sander_output.py',
##               'drivers/run_ptraj.py',
##               'drivers/write_ptraj_input_files.py',
##               'drivers/setup_hbond_ptraj.py',
##               'drivers/combine_hbonds.py',
##               'drivers/compare_hbonds.py',
##               'drivers/subset_hbonds.py',
##               ],
      )
