Correlated Dynamics
===================

Executive Summary
-----------------

Here's how to make plots for a simulation of 1rx1 with 1000ps windows:

Make the directories
~~~~~~~~~~~~~~~~~~~~

::

  ssh node6
  cd /data/people/mlerner
  mkdir correlated_dynamics
  cd correlated_dynamics
  mkdir ptraj_files
  mkdir ptraj_files/images
  mkdir ptraj_files/1rx1

Calculating the correlation matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write_ptraj_input_files.py --strip-water --strip-hydros --input-dir=. --mdcrd=1rx1.trj --structure-name=1rx1 --start=100 --stop=1000 --window-size=200 --window-spacing=100

::

  #
  # Standard defaults are windowsize of 1000ps (1NS) and windowspacing of 100ps.  
  #
  write_ptraj_input_files.py --strip-water --strip-hydros --input-dir=. --mdcrd=1rx1.mdcrd 
     --structure-name=1rx1 --start=500 --stop=5500 --ps=5 --window-size=1000 --window-spacing=500
  
  nohup run_ptraj.py --input-dir=. --prmtop=1rx1.prmtop --structure=1rx1&
  
  cd ptraj_files/1rx1
  bzip2 *.dat #*
  cd ../..

Extract the per-residue and per-atom information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  nohup do_correlated_md_analysis.py --output-dir=ptraj_files --structure-name=1rx1 --start=500 
     --stop=5500  --window-size=1000 --non-ca-resis=160,161-175&

Make the plots
~~~~~~~~~~~~~~

::

  nohup make_correlated_dynamics_plots.py --structure-name=1rx1 --start=500 --stop=5500 
     --window-size=1000 --window-spacing=500&
  
  make_movies.py --structure-name=1rx1 --plot-types=ca,avg,max,min,abs,straight,mainheavy
  
  cd ptraj_files
  tar cvf 1rx1_100ps_movies.tar 1rx1BigAnimatedMovies.html images/animated_*
  mv 1rx1_100ps_movies.tar ~
  cd ..


More detailed explanations
--------------------------

The scripts have many options, and we'll describe a standard setup
here. You'll need a couple of things before you begin:

 1. An MD trajectory from `sander`. In this example, it will be called
    `1rx1.mdcrd`.

 2. The paramater/topology file corresponding to that trajectory. Ours
    will be called `1rx1.prmtop`

 3. Sufficient disk space and processor power. This can easily eat up
    several gigs of disk space, and we usually run things with 2G of
    memory.

Making the directories
~~~~~~~~~~~~~~~~~~~~~~

First of all, we need to set up a directory structure to store our
files. We need a main directory which will contain all of our results
(`correlated_dynamics`). Inside that, we need to store the images and
the data. We'll store images and data in `ptraj_files`. Images will go
in `ptraj_files/images`. We'll have different directories for each
different structure that we study. These examples will be for the DHFR
structure 1RX1, and we'll store data in `ptraj_files/1rx1`. Assuming
we do all of this on node6 of our cluster, here's how to set up the
directories::


  ssh node6
  cd /data/people/mlerner
  mkdir correlated_dynamics
  cd correlated_dynamics
  mkdir ptraj_files
  mkdir ptraj_files/images
  mkdir ptraj_files/1rx1


Calculating the correlation matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We use `ptraj` to calculate the correlation matrices. So, we need to
write out lots of ptraj input files and then we need to run them. This
sets up all of the calculations that we'll use, so it's important to
get the options right.


 * `structure-name` should be the same thing you specified
   earlier. `1rx1` in our case.

 * `input-dir` is the path to the directory that contains the
   trajectory and parameter/topology file.

 * `strip-water` and `strip-hydros` tell us whether or not to include
   waters and hydrogens in our calculations. Waters are almost never
   worth including. Hydrogens can be interesting, but it's very easy
   to run out of memory on reasonably-sized systems, so we typically
   exclude them.

 * `ps` tells us how often (in picoseconds) frames were written to the
   mdcrd file.

 * `start` and `stop` tell us (in ps) when to start and stop the
   windows. `window-size` tells us how long each window should be (in
   ps). `window-spacing` tells us how often to write out windows (in
   ps). The defaults are to write out windows of length 1ns (1000ps)
   every 100ps.

 * There are several other options, including the ability to insert
   user-specified commands directly to the ptraj input files. Please
   use the `help` option for more information.

Finally, these files can take up a lot of disk space. We typically
compress things with bzip2. The scripts are smart enough to decompress
things on the fly later on.

so, here's how we set up and run ptraj::


  #
  # Standard defaults are windowsize of 1000ps (1NS) and windowspacing of 100ps.  
  #
  write_ptraj_input_files.py --strip-water --strip-hydros --input-dir=. --mdcrd=1rx1.mdcrd 
     --structure-name=1rx1 --start=500 --stop=5500 --ps=5 --window-size=1000 --window-spacing=100
  
  nohup run_ptraj.py --input-dir=. --prmtop=1rx1.prmtop --structure=1rx1&
  
  cd ptraj_files/1rx1
  bzip2 *.dat
  cd ../..


Extract the per-residue and per-atom information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ptraj calculates correlations between each atom. We also want to
calculate the following quantities:

 * Correlations between alpha-carbons.

 * Correlations between main-chain heavy atoms

 * The following quantities on a per-residue basis:

   * average

   * maximum

   * minimum

   * largest absolute value

If hydrogens are included, we will also calculate correlations between
potential hydrogen-bond donors and acceptors.

Since we are calculating alpha-carbon correlations, it is important to
provide a list of residues that do not contain alpha carbons
(`non-ca-resis`). In this case, 160 is our cofactor and 161-175 are
our counter-ions. We could have excluded these during the
`write_ptraj_input_files.py` command, but chose not to.

All of our files follow a standard naming convention, so telling each
successive command the start, stop, spacing and size of the windows is
enough to make sure that the correct files are read in.

::

  nohup do_correlated_md_analysis.py --input-dir=. --structure-name=1rx1 --start=500 --stop=5500  
     --window-size=1000 --non-ca-resis=160,161-175&


Make the plots
~~~~~~~~~~~~~~

Now we have calculated everything and it's time to make the
plots. `make_correlated_dynamics_plots.py` makes the individual plots
and has *many* different options (again, use `help` to list them
all). Here is a standard run.


 * we specify the structure and the details about the windows as
   before.

 * `plot-types` is a comma-separated list of plot types. The standard
   ones that we use are `ca,avg,max,min,abs,straight,mainheavy` and
   the command-line help will detail other options for you. Since we
   don't specify this on the command-line below, it will default to
   the standard options.

::

  nohup make_correlated_dynamics_plots.py --input-dir=. --structure-name=1rx1 --start=500 
     --stop=5500 --window-size=1000&


That produces plots of each of the individual windows. It's worth
examining these on their own. However, it's often a lot more
interesting to look at movies of these all pasted together. The
`convert` program is used to do this. If it's not installed, it's easy
to install on OS X, Linux and Windows. `make_movies.py` calls
`convert` appropriately::

  make_movies.py --structure-name=1rx1 --plot-types=ca,avg,max,min,abs,straight,mainheavy


Finally, we may wish to move the movies to another machine for
viewing. The movies are in the `images` subdirectory. `make_movies.py`
also generates an html file that shows all of the movies with
thumbnails. Here's how to collect the movies and html file::

  cd ptraj_files
  tar cvf 1rx1_100ps_movies.tar 1rx1BigAnimatedMovies.html images/animated_*
  mv 1rx1_100ps_movies.tar ~
  cd ..


Other notes
~~~~~~~~~~~

Specific options
++++++++++++++++

`make_correlated_dynamics_plots.py` has several useful options.

 * Residues of interest can be marked with `mark-resis`. This is
   useful in several cases, including:

   * marking loops, helices or other regions of interest

   * marking every 10th residue to show a grid (especially useful for
     orientation in the main-chain heavy plots)

   * you can select different color maps with the `cmap` option.

Please read through the command-line help for a detailed, up-to-date explanation.

File formats and how to save space
++++++++++++++++++++++++++++++++++

The bzip2'd files are very small. However, especially when dealing
with main-chain heavy atoms, they can be quite slow. `numpy` has an
internal format that is much faster to read. You can use
`convert_to_numpy_format.py` to covert your output to this format. It
looses a small amount of precision, but that seems to be completely
negligible. All of the correlated-dynamics scripts will deal
transparently with any combination of .dat, .numpy and .bz2 files. You
specify things as above and it'll figure out how to read .dat or
.numpy or .numpy.bz2 without any trouble.

Smaller movies
++++++++++++++

The movie files are created as animated GIFs. This has the advantage
that they can be played anywhere. However, they can be quite
large. Sophisticated tools such as ``mencoder``
(http://www.mplayerhq.hu) can use two-pass encoding to convert them
into much smaller AVI files. We find that the XVid codec produces good
quality movies. Codecs and encoders are frequently updated, and it is
suggested that users consult websites like the one mentioned above for
the most current information.

Documentation from the scripts
------------------------------

write_ptraj_input_files.py
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  Usage: write_ptraj_input_files.py [options]
  
  Please make sure that you have created the following directories:
  
    output-dir
    output-dir/structure-name
    output-dir/images
  
  This script will emit the PDB file that you will use as a
  reference structure later on.  It will live in <outputdir>/<structure>_ref.pdb.1
  
  
  Options:
    -h, --help            show this help message and exit
    --structure-name=STRUCTURENAME
                          Name of your structure.  E.g. 1RX1 or 1SGZ. [default:
                          1rx1]
    --output-dir=OUTPUTDIR
                          Directory where we will put our results.  This should
                          be the same as the directory where we put our ptraj
                          files before, and it should contain the .dat files
                          that ptraj outputs.  [default: ptraj_files/]  The
                          images will go to <outputdir>/images and the html file
                          will be in <outputdir>
    --start=START         Time, in ps, to start the windows.  [default: 500]
    --stop=STOP           Time, in ps, to stop the windows.  [default: 10500]
    --window-size=WINDOWSIZE
                          Length, in ps, of window size.  [default: 1000]
    --window-spacing=WINDOWSPACING
                          Spacing between windows, in ps.  [default: 100]
    --input-dir=INPUTDIR  Directory that contains our input files.  It should
                          contain the prmtop file and the mdcrd file.  [default:
                          ./]
    --mdcrd=MDCRD         Comma-separated list of mdcrd files
    --ps=PS               Number of ps per frame.  [default: 5]
    --align=ALIGN         How to align. 'all' means 'rms first *'.  'none' means
                          no alignment.  Any other string will be treated as the
                          alignment string.  For example, if you say ':1-428@CA'
                          the ptraj file will say 'rms first :1-428@CA'.
                          [default: all]
    --strip-hydros        Strip the hydrogens out during the ptraj runs.
    --strip-waters        Strip the waters out during the ptraj runs.  This
                          assumes they're named WAT.
    --write-covar         Write out the covariance matrix [default: False]
    --other-ptraj-strips=OTHER_PTRAJ_STRIPS
                          Comma-separated list of other things that ptraj should
                          strip.  For example, could be :WAT,:BOB and we would
                          add two lines, one saying strip :WAT and one saying
                          strip :BOB.

run_ptraj.py
~~~~~~~~~~~~

::

  Usage: run_ptraj.py [options]
  
  Please make sure that you have created the following directories:
  
    output-dir
    output-dir/structure-name
    output-dir/images
  
  
  Options:
    -h, --help            show this help message and exit
    --structure-name=STRUCTURENAME
                          Name of your structure.  E.g. 1RX1 or 1SGZ. [default:
                          1rx1]
    --output-dir=OUTPUTDIR
                          Directory where we will put our results.  This should
                          be the same as the directory where we put our ptraj
                          files before, and it should contain the .dat files
                          that ptraj outputs.  [default: ptraj_files/]  The
                          images will go to <outputdir>/images and the html file
                          will be in <outputdir>
    -i INPUTDIR, --input-dir=INPUTDIR
                          Directory that contains our input files.  It should
                          contain the prmtop file and the mdcrd file.  [default:
                          ./]
    -p PRMTOP, --prmtop=PRMTOP
                          Name of prmtop file
  
do_correlated_md_analysis.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  Usage: do_correlated_md_analysis.py [options]
  
  Please make sure that you have created the following directories:
  
    output-dir
    output-dir/structure-name
    output-dir/images
  
  
  Options:
    -h, --help            show this help message and exit
    --structure-name=STRUCTURENAME
                          Name of your structure.  E.g. 1RX1 or 1SGZ. [default:
                          1rx1]
    --output-dir=OUTPUTDIR
                          Directory where we will put our results.  This should
                          be the same as the directory where we put our ptraj
                          files before, and it should contain the .dat files
                          that ptraj outputs.  [default: ptraj_files/]  The
                          images will go to <outputdir>/images and the html file
                          will be in <outputdir>
    --start=START         Time, in ps, to start the windows.  [default: 500]
    --stop=STOP           Time, in ps, to stop the windows.  [default: 10500]
    --window-size=WINDOWSIZE
                          Length, in ps, of window size.  [default: 1000]
    --window-spacing=WINDOWSPACING
                          Spacing between windows, in ps.  [default: 100]
    --non-ca-resis=NON_CA_RESIS
                          Comma separated list of residues that don't contain
                          alpha carbons.  We need this to make some of our
                          output images.  [default: []], but you could say
                          160,161 for example
  
make_correlated_dynamics_plots.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  Usage: make_correlated_dynamics_plots.py [options]
  
  Please make sure that you have created the following directories:
  
    output-dir
    output-dir/structure-name
    output-dir/images
  
  This will spit out an html file that will show you your images.  If you
  need to look at the images on another machine, tar up the html file and
  output-dir/images together and move that to the other machine.
  
  
  Options:
    -h, --help            show this help message and exit
    --structure-name=STRUCTURENAME
                          Name of your structure.  E.g. 1RX1 or 1SGZ. [default:
                          1rx1]
    --output-dir=OUTPUTDIR
                          Directory where we will put our results.  This should
                          be the same as the directory where we put our ptraj
                          files before, and it should contain the .dat files
                          that ptraj outputs.  [default: ptraj_files/]  The
                          images will go to <outputdir>/images and the html file
                          will be in <outputdir>
    --start=START         Time, in ps, to start the windows.  [default: 500]
    --stop=STOP           Time, in ps, to stop the windows.  [default: 10500]
    --window-size=WINDOWSIZE
                          Length, in ps, of window size.  [default: 1000]
    --window-spacing=WINDOWSPACING
                          Spacing between windows, in ps.  [default: 100]
    --cmap=CMAP           Color map to use when making the plots. Our custom
                          cmaps are Normal and Scaled.  Standard matplotlib
                          cmaps ['Spectral', 'summer', 'RdBu', 'gist_earth',
                          'Set1', 'Set2', 'Set3', 'Dark2', 'hot', 'RdPu',
                          'YlGnBu', 'RdYlBu', 'gist_stern', 'cool', 'gray',
                          'GnBu', 'gist_ncar', 'gist_rainbow', 'bone', 'RdYlGn',
                          'spring', 'Accent', 'PuBu', 'spectral', 'gist_yarg',
                          'BuGn', 'YlOrRd', 'Greens', 'PRGn', 'gist_heat',
                          'Paired', 'hsv', 'Pastel2', 'Pastel1', 'copper',
                          'OrRd', 'jet', 'BuPu', 'Oranges', 'PiYG', 'YlGn',
                          'gist_gray', 'flag', 'BrBG', 'Reds', 'RdGy', 'PuRd',
                          'Blues', 'Greys', 'autumn', 'pink', 'binary',
                          'winter', 'prism', 'YlOrBr', 'Purples', 'PuOr',
                          'PuBuGn'] are also supported.[default: Normal]
    --plot-types=PLOTTYPES
                          Comma-separated list of plot types.  [default: ['ca',
                          'avg', 'max', 'min', 'abs', 'straight', 'mainheavy',
                          'allheavy', 'sidechainhbond', 'hbond']]
    --mark-resis=MARKRESIS
                          A list of residues to mark on the plots.  [default:
                          []]
    --highlight=HIGHLIGHT
                          How strongly to highlight the marked residues.  Note
                          that --highlight-mode tells us how exactly we will do
                          the highlighting.  0.1 and 0.2 are decent values if
                          you want to use this feature for most plots, although
                          you'll need something stronger for the absolute value
                          plots. [default: 0.2]
    --highlight-mode=HIGHLIGHTMODE
                          When highlight-mode is 'negative' we put a white block
                          down on top of the marked residues, the opacity of
                          which is controled by --highlight.  When it's
                          'positive', we put that white block down on squares of
                          residues that are *not* highlighted instead.  When
                          it's 'supernegative', we do just like 'negative'
                          except that the block will be twice as opaque where
                          the highlighted rows and columns intersect.  Positive
                          and supernegative seem to be more useful than
                          negative.  [default: positive]
    --skip-resis=SKIPRESIS
                          A list of residues that will be skipped in the plots.
                          [default: []]
    --no-ticks            do not include tick marks on the axes
    --dpi=DPI             dpi for figures [default: 200]
    --title=TITLE         If no title is specified, one will be automatically
                          generated. Note that the title is part of the filename
                          that we save. [default: none]

make_movies.py
~~~~~~~~~~~~~~

::

  Usage: make_movies.py [options]
  
  Please make sure that you have created the following directories:
  
    output-dir
    output-dir/structure-name
    output-dir/images
  
  
  Options:
    -h, --help            show this help message and exit
    --structure-name=STRUCTURENAME
                          Name of your structure.  E.g. 1RX1 or 1SGZ. [default:
                          1rx1]
    --output-dir=OUTPUTDIR
                          Directory where we will put our results.  This should
                          be the same as the directory where we put our ptraj
                          files before, and it should contain the .dat files
                          that ptraj outputs.  [default: ptraj_files/]  The
                          images will go to <outputdir>/images and the html file
                          will be in <outputdir>
    --plot-types=PLOTTYPES
                          Comma-separated list of plot types.  [default: ['ca',
                          'avg', 'max', 'min', 'abs', 'straight', 'mainheavy',
                          'allheavy', 'sidechainhbond', 'hbond']]
    --no-slow-movies      Set this if you do not want to generate the movies
                          that have 0.5s spacing between the frames.
    --movie-link=MOVIELINK
                          'fast' if you want the thumbnails to link to the fast
                          images, anything else for the slow ones. [default:
                          fast]
  
convert_to_numpy_format.py
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  Usage: This will convert the .dat or .dat.bz2 files to numpy versions.
      It will not automatically delete the .dat(.bz2) files.  If your input
      files are bz2, your output files will be too.
  
      If you already have a corresponding .numpy or .numpy.bz2 file, we won't
      write out a new file.
      
  
  Options:
    -h, --help            show this help message and exit
    --dir=DIR             Directory in which the files reside. [default: .]
    --structure-name=STRUCTURENAME
                          Name of your structure.  E.g. 1RX1 or 1SGZ
    --compression=COMPRESSION
                          Type of compression currently used on files.  Leave
                          blank for uncompressed, .bz2 if they're .bz2 files.
                          Note that it's '.bz2' not 'bz2'. [default: ]
    --all-dat-files       By default, we will only convert the
                          all_atom_correlmat files.  If you use this option, we
                          will convert all dat files.  Don't forget, though,
                          that you'll still have to call this command twice if
                          you have some files that are bzipped and some that are
                          not.
  
