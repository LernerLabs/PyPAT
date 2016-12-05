Hydrogen Bonding
================

Executive Summary
-----------------

A commonly performed analysis of molecular dynamics trajectories is
hydrogen bond analysis, in which statistics of all the hydrogen bonds
in a system (occupancy, lifetime, distance, etc.) are monitored
throughout the simulation.  However, this type of analysis is
difficult for two reasons: it requires very large amounts of computer
memory, and the abundance of data obtained can be overwhelming.  The
PyPAT H-bond analysis tools provide solutions to these issues.

The tools come as a suite of four programs: ``setup_hbond_ptraj.py``,
``combine_hbonds.py``, ``subset_hbonds.py``, and
``compare_hbonds.py``. The first two programs provide a workaround to
the memory issue, and the last two provide an easy way to comb through
and sort the data.

The programs make use of the trajectory analysis program of the AMBER
software package, ptraj.  Because of memory requirements, ptraj can
only be used to perform a global hydrogen-bonding analysis on a
reasonably-sized system for short lengths (several hundred ps) of a
trajectory.  In order to analyze a long trajectory, we first divide it
into short segments, use ptraj to analyze each one, and then recompile
the data into a unified set.  These tasks are accomplished by
``setup_hbond_ptraj.py`` and ``combine_hbonds.py``.  The other two
programs provide methods to select a desired subset of data, sort it,
and output it in a clean-looking format.

Typically, a series of commands similar to the following will be used
in the analysis:

::

  $ setup_hbond_ptraj.py  -x 1sgz.mdcrd -p 1sgz.prmtop -b 1 -e 5000 -g 500 -n 389  

  $ ./run_hbond_ptraj  

  $ combine_hbonds.py segment*.out -p 1sgz.prmtop -r 4 -o combined1.out  

  $ subset_hbonds.py combined1.out -R 2-6,9-12 -O 60 -y -c   

  $ compare_hbonds.py combined1.out combined2.out -i 1sgz,1w50 -R 9-27 -A bb_only -O 20 -s donor -y -c 

The operation of each of these programs is detailed below. We note
that most of these scripts take the ``-D`` option, which tells the
scripts that the data files live in a different directory than the
current one.

setup_hbond_ptraj.py
--------------------

This program prepares a series of input scripts that can be run by
ptraj as well as a script that will automatically perform the
execution of ptraj with each one.

Documentation from the script::

  Usage: setup_hbond_ptraj.py [options]
  This program sets up ptraj input files to perform a series of H-bond 
  analyses on a trajectory.  The trajectory is broken into segments of a 
  specified size.  Files written out include a ptraj input file for each 
  segment and a file "run_hbond_ptraj" that will sequentially execute
  ptraj with each one.  
  
  The only required option is --mdcrd-file (-x), which specifies the
  coordinate file.  
  
  Important notes:  
  * A file called "mask" can be created to include hydrogen-bonding
    residues other than the standard amino acids and water.  The lines in
    mask will be copied "as is" into each ptraj input file, so follow the
    format for specifying hydrogen bond donors and acceptors in the ptraj
    documentation.
  
  Options:
    -h, --help            show this help message and exit
    -x MDCRD_FILE, --mdcrd-file=MDCRD_FILE
                          Coordinate file to use for H-bond analysis (REQUIRED).
                          [default: none]
    -p PRMTOP_FILE, --prmtop-file=PRMTOP_FILE
                          Amber parameter/topology file to use for H-bond
                          analysis.  If None, the name will be guessed by
                          replacing the extension of the coordinate file name
                          with ".prmtop"  [default: none]
    -D MDCRD_PRMTOP_DIR, --mdcrd-prmtop-dir=MDCRD_PRMTOP_DIR
                          The directory that contains the coordinate and prmtop
                          file.  If None, the names of the files will be used
                          without modification.  [default: none]
    -o OUTPUT_FILE_BASE, --output-file-base=OUTPUT_FILE_BASE
                          The first part of the name of the ptraj input files
                          that will be created.  The files will be named
                          OUTPUT_FILE_BASEnn.in, where nn is the two-digit
                          segment number.  Once ptraj is run (with
                          run_hbond_ptraj), the output files containing the
                          H-bond data will be named OUTPUT_FILE_BASEnn.out.
                          [default: segment]
    -B HBOND_DATA_DIR, --hbond-data-dir=HBOND_DATA_DIR
                          The directory where the ptraj input files will be
                          placed.  [default: .]
    -b BEGIN_FRAME, --begin-frame=BEGIN_FRAME
                          The number of the first frame to use in the H-bond
                          analysis.  [default: 1]
    -e END_FRAME, --end-frame=END_FRAME
                          The number of the last frame to use in the H-bond
                          analysis.  [default: 10000]
    -g SEGMENT_SIZE, --segment-size=SEGMENT_SIZE
                          The number of frames included in each segment of the
                          trajectory.  [default: 1000]
    -n NUM_RESI, --num-resi=NUM_RESI
                          The number of residues in the protein.  [default:
                          10000]
    -d DIST_CUTOFF, --dist-cutoff=DIST_CUTOFF
                          The distance cutoff that defines whether or not an
                          interaction is an H-bond.  [default: 3.0]
    -a ANGLE_CUTOFF, --angle-cutoff=ANGLE_CUTOFF
                          The angle cutoff that defines whether or not an
                          interaction is an H-bond.  [default: 120.0]
    -s, --no-self         Flag to turn off the inclusion of H-bonds between
                          atoms within the same residue.  [default: True]
    -S, --solvent         Flag to include solvent-protein interactions in the
                          analysis.  [default: False]
    -m MASK_FILE, --mask-file=MASK_FILE
                          File that contains extra lines to include in the ptraj
                          input files, primarily to include masks for ligands.
                          [default: mask]
  

The only required input to ``setup_hbond_ptraj.py`` is the coordinate
file given with the ``-x`` (``--mdcrd-file``) option.  The program
will terminate with an error if this input is not given.  In addition,
the program will terminate with an error if the specified coordinate
or prmtop file does not exist.

The program creates a number of files in the directory specified with
the ``-B`` (``--hbond-data-dir``) option.  The trajectory is broken
into several segments, beginning with the frame specified with the
``-b`` (``--begin-frame``) option, ending with the frame specified
with the ``-e`` (``--end-frame``) options, and with each segment
containing a number of frames specified by the ``-g``
(``--segment-size``) option.  The number of segments (and the number
of ptraj input files created) is::

  (END_FRAME - BEGIN_FRAME + 1) / SEGMENT_SIZE.

For statistical purposes, the segments must all have the same size, so
if ``SEGMENT_SIZE`` initially does not divide the numerator equally,
``END_FRAME`` is decreased so that it does.  A warning will be
presented to the user to indicate that frames are being removed from
consideration.  The different ptraj input files that are created are
identical except for the particular frames of the trajectory that they
will be used to analyze.  These files are named based on the ``-o``
(``--output-file-base``) option; the names are
``OUTPUT_FILE_BASEnn.in``, where nn is the two-digit segment number.

An additional file called ``run_hbond_ptraj`` is also created.  This
is an executable file that should be run immediately after
``setup_hbond_ptraj.py`` is finished.

The program can handle all of the amino acid residues recognized by
the Amber programs tLEaP and xLEaP (one of which, presumably, was used
to prepare the system for MD).  These include the twenty natural amino
acids and the following additional residues: His in its delta-,
epsilon-, and doubly-protonated forms (named HID, HIE, and HIP,
respectively); neutralized Asp, Glu, and Lys (named ASH, GLH, and
LYN); and Cys in its disulfide and deprotonated forms (named CYX and
CYM).  If the system of interest contains residues that are not among
these amino acids, such as ligands, the additional hydrogen-bond
donors and acceptors can be specified in a file named mask [or an
alternate name specified by the ``-m`` (``--mask-file``) option].  The
mask file will be read line for line into each ptraj input file, so it
must contain the appropriate syntax for specifying the donors and
acceptors to ptraj.  In ptraj, hydrogen bond "donors" are defined as
the heavy atoms that are not covalently bound to the hydrogen atom,
and the "acceptors" are the heavy atoms that are covalently bound to
the hydrogen.  This is the opposite definition of the normal usage of
the words, but it is the convention that ptraj has adopted.  The file
mask should contain one line for each potential H-bond donor and
acceptor according to the syntax::

  donor mask :lig@atom-name
  acceptor mask :lig@heavy-atom-name :lig@H-atom-name

where ``lig`` is the ligand residue name or number and the atom-names
are the names of the participating atoms.

The ``-n`` (``--num-resi``) option to specify the number of residues
is of no consequence if there is no solvent present in the system of
interest, but it is important if there is solvent.  With solvent
molecules present, failing to provide the correct number of residues
will result in some of the water oxygen atoms being treated explicitly
as hydrogen bond acceptors.  This will not affect the analysis of
protein-protein hydrogen bonds, but it will result in the compilation
of more data than is necessary and may result in the memory issues
during the ptraj execution.  If the ``-n`` option is not given, a
warning will be printed.  In some cases, the inclusion of H-bonds
between protein and solvent may be desirable, but the ``-S``
(``--solvent``) flag should be used instead of increasing the number
of residues.  See the ptraj documentation for details on the way
solvent donors and solvent acceptors are handled.

The ``-d`` (``--dist-cutoff``), ``-a`` (``--angle-cutoff``), ``-s``
(``--no-self``), and ``-S`` (``--solvent``) options can be used to
control what interactions ptraj considers to be hydrogen bonds.  The
default distance and angle cutoffs are the same as those of ptraj, but
in contrast to ptraj, hydrogen bonds between atoms of the same residue
*will* be reported unless this behavior is turned off with the ``-s``
flag.

run_hbond_ptraj
---------------

This program is created by ``setup_hbond_ptraj.py`` to perform the
ptraj executions.  It takes no arguments.  ``run_hbond_ptraj``
sequentially performs the ptraj execution for each segment.  The
resulting files output by ptraj have the name
``OUTPUT_FILE_BASEnn.out``, again where ``nn`` is the two-digit
segment number.

The output files contain H-bond data.  For each H-bond, ptraj reports
the occupancy percentage, average heavy atom-heavy atom distance,
average heavy atom-H-heavy atom angle, average lifetime, and the
maximum number of continuous frames the H-bond is populated.  Standard
deviations of the distance, angle, and lifetime are also included.  In
addition, a 10-character "graph" that displays the occupancy in each
tenth of the trajectory is reported.  Each character indicates a
different level of occupancy: a space (0-5%), . (5-20%), - (20-40%), o
(40-60%), x (60-80%), * (80-95%), or @ (95-100%).

combine_hbonds.py
-----------------

Once the files with the H-bond data have been generated by ptraj,
combine_hbonds.py is used to compile the data into a single file.

Documentaiton from the script::

  Usage: combine_hbonds.py FILE1 [ FILE2 [ ... ] ] [options] 
  FILE1 and additional optional FILEs are files containing hbond data that
  were produced by the ptraj hbond command.  At least one such file is
  required.  The data in the files are spliced together to created a 
  unified data set.  
  
  Important notes:
  * An AMBER prmtop or a PDB file may be input with the -p option.  The
    file will be used to determine the residue name associated with each
    residue number in the system.  If the file name does not end with
    '.pdb', the file will be assumed to be an AMBER prmtop file.  The 
    offset and amino acid code are also used in the residue name 
    generation.
  * The resi_criteria option takes a comma-separated string containing any
    or all of the following:
          - individual residue numbers
          - a range of numbers, separated by a '-'
          - strings associated with valid residue lists in the
              standard file 'residue_lists' 
  * The atom_criteria option takes a comma-separated string containing
    the atom names to report.  In addition, the string can contain any
    of these strings:
          - 'bb_only': only H-bonds between two backbone atoms
          - 'not_bb': no H-bonds between two backbone atoms
          - 'protein_only': no H-bonds involving water
  
  Options:
    -h, --help            show this help message and exit
    -o OUTPUT_FILE, --output-file=OUTPUT_FILE
                          The name of the output file.  If None, the results
                          will be written to stdout.  Supplying a name is
                          recommended.  [default: none]
    -B HBOND_DATA_DIR, --hbond-data-dir=HBOND_DATA_DIR
                          The directory that contains the H-bond data files.  If
                          None, the file names will be used without
                          modification, and output will be written to the
                          current directory.  [default: none]
    -g SEGMENT_SIZE, --segment-size=SEGMENT_SIZE
                          The number of frames included in each segment of the
                          trajectory.  [default: 1000]
    -p PRMTOP_FILE, --prmtop-file=PRMTOP_FILE
                          Amber parameter/topology file or PDB file for the
                          system.  [default: none]
    -D PRMTOP_DIR, --prmtop-dir=PRMTOP_DIR
                          The directory that contains the prmtop (or PDB) file.
                          If None, the name of the file will be used without
                          modification.  [default: none]
    -r RESI_OFFSET, --resi-offset=RESI_OFFSET
                          The offset between the residue numbers in the prmtop
                          (or PDB) file and the actual residue numbers.
                          [default: 0]
    -a AA_CODE, --aa-code=AA_CODE
                          Must be 1 or 3.  Indicates the use of 1- or 3-letter
                          amino acid codes in residue names.  [default: 3]
    -y, --occ-graph-only  Flag to report only the occupancy and graph data (no
                          distance or angle data).  [default: False]
    -R RESI_CRITERIA, --resi-criteria=RESI_CRITERIA
                          A comma- and dash-separated list of residue numbers to
                          include in the analysis.  [default: all]
    -A ATOM_CRITERIA, --atom-criteria=ATOM_CRITERIA
                          A comma-separated list of atom names to include in the
                          analysis.  [default: all]
    -O OCC_THRESH, --occ-thresh=OCC_THRESH
                          The minimum occupancy threshold that the H-bonds must
                          have to be reported.  [default: 0.0]


``FILE1`` and additional optional ``FILE`` s are the files containing
H-bond data that were produced by the ptraj hbond command.  At least
one such file is required.  If the data files are not contained in the
current directory, the ``-B`` (``--hbond-data-dir``) option can be
used to indicate where they are located.  The name of the output file
can be specified by the ``-o`` (``--output-file``) option.  Because a
file output from ``combine_hbonds.py`` is required for the analysis
aids subset_``hbonds.py`` and ``compare_hbonds.py`` (below), it is
recommended that a name be supplied.  The output file will also be
located in directory ``HBOND_DATA_DIR``.

A line is output providing information for each H-bond in the system.
As shown in the diagram below, this information includes the
participating atoms, the occupancy percentage, the total number of
frames included in the analysis, the average distance and angle of the
H-bond (and their standard deviations), and the graph.  Lifetime and
maximum continuous occupancy data are not included, because splitting
the trajectory into segments artificially shortens these values.
Consequently, they are not reliable indicators of the data over the
trajectory as a whole.

::

  Thr376  OG1--HG1 ... OG   Ser295   35.20(  250) 2.854(0.17) 22.49(12.08) |*x@.    x@|@@.--.  . |
  ----- participating atoms ------   -occ- -num-  -dist -SD-  angle --SD-   ------- graph -------
                                      pct  frames


Several options control the display of the output.  An AMBER prmtop or
a PDB file may be input with the ``-p`` (``--prmtop-file``) option.
(If the file name does not have an extension of ``.pdb``, the program
will assume it is a prmtop file.)  The file will be used to determine
the residue type of each residue number in the system.  The residue
type is prepended to the residue number to create the full residue
name displayed in the output file.  Whether to use a three- or
one-letter amino acid code is determined by the input to the ``-a``
(``--aa-code``) option.  Additionally, an offset can be supplied to
change the numbering of the residues using the "r" (--resi-offset)
option.  This may be useful in systems with an unusual numbering
convention.  For example, the residues of the 1SGZ crystal structure
are numbered from -3 to 385, but the Amber setup shifts them to 1-389.
Using an offset of -4 brings the numbering back in line with
convention.  If a prmtop file for 1SGZ and this offset are given to
the program, the residues in the output are labeled, for instance,
"Tyr71" instead of "75."  Also, the ``-y`` (``--occ-graph-only``) flag
will affect the output presentation of the H-bond data.  Using this
flag will prevent the reporting of the distance and angle data, so
that only the occupancy percentage and graph will be displayed.

The data in the input files can be filtered prior to output by
specifying residue or atom criteria or an occupancy threshold.  The
residue criteria (``-R``, ``--resi-criteria``) option takes a
comma-separated string containing the residues numbers to report.  It
may contain any or all of the following: individual residue numbers, a
range of numbers separated by a dash, or strings associated with valid
residue lists in the standard file residue_lists as described below.
The atom criteria (``-A``, ``--atom-criteria``) option takes a
comma-separated string containing the atom names to report.  In
addition, the program understands these strings: "``bb_only``", which
will include only H-bonds between two backbone atoms; "``not_bb``",
which will exclude H-bonds between two backbone atoms; or
"``protein_only``", which will exclude any H-bonds involving water.
In addition, both the residue and atom criteria can take the word
"``all``".  Lastly, the occupancy threshold (``-O``, ``--occ-thresh``)
option will exclude all H-bonds with an occupancy percentage below the
threshold.

Note on the ``residue_lists`` file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The lines of the ``residue_lists file`` can be used to associate a
list of residues with a particular name according to the following
syntax:

::

  name:  residue_string

where ``name`` is any string and ``residue_string`` is any comma- and
dash-separated list of residues.  For example, the line

::

  loops: 9-12,65-67

will associate the string ``loops`` with the residue list 9, 10, 11,
12, 65, 66, and 67.  A residue list is even allowed to include strings
that were defined above it.  For example, if the following lines are
contained in the residue_lists file:

::

  loops     : 9-12,65-67
  more_loops: loops,100,104

the name ``more_loops`` is associated with residues 9, 10, 11, 12, 65,
66, 67, 100, and 104.

subset_hbonds.py
----------------

Once a file of combined H-bond data is created with
``combine_hbonds.py``, ``subset_hbonds.py`` allows the user to select
and sort a subset of the data, which greatly facilitates analysis.

Documentation from the script::

  Usage: subset_hbonds.py FILE1 [options] 
  FILE1 is a file created by combine_hbonds.py that contains a dataset of
  H-bonds.  Only a subset of all the data is presented according to the
  criteria presented by the user.  The data will also be sorted according
  to the metric specified by the user. 
  
  Important notes:
  * The resi_criteria option takes a comma-separated string containing any
    or all of the following:
          - individual residue numbers
          - a range of numbers, separated by a '-'
          - strings associated with valid residue lists in the
              standard file 'residue_lists'
  * The atom_criteria option takes a comma-separated string containing
    the atom names to report.  In addition, the string can contain any
    of these strings:
          - 'bb_only': only H-bonds between two backbone atoms
          - 'not_bb': no H-bonds between two backbone atoms
          - 'protein_only': no H-bonds involving water
  * Note that ptraj uses a definition of H-bond donor and acceptor that is
    opposite of the normal convention.  This program follows the
    definitions of ptraj, in which the acceptor is the atom covalently 
    bonded to the hydrogen atom.
  
  Options:
    -h, --help            show this help message and exit
    -o OUTPUT_FILE, --output-file=OUTPUT_FILE
                          The name of the output file.  If None, the results
                          will be written to stdout.  [default: none]
    -B HBOND_DATA_DIR, --hbond-data-dir=HBOND_DATA_DIR
                          The directory that contains the H-bond data files.  If
                          None, the file names will be used without
                          modification, and output will be written to the
                          current directory.  [default: none]
    -s SORT, --sort=SORT  The quantity used to sort the results.  Must be one of
                          "occ_pct" (occupancy percentage), "donor", "acceptor",
                          "dist", or "angle".   [default: occ_pct]
    -c, --compress        Flag to compress the H-bond graph.  [default: False]
    -y, --occ-graph-only  Flag to report only the occupancy and graph data (no
                          distance or angle data).  [default: False]
    -R RESI_CRITERIA, --resi-criteria=RESI_CRITERIA
                          A comma- and dash-separated list of residue numbers to
                          include in the analysis.  [default: all]
    -A ATOM_CRITERIA, --atom-criteria=ATOM_CRITERIA
                          A comma-separated list of atom names to include in the
                          analysis.  [default: all]
    -O OCC_THRESH, --occ-thresh=OCC_THRESH
                          The minimum occupancy threshold that the H-bonds must
                          have to be reported.  [default: 0.0]

The only requirement for this program is a single argument ``FILE1``,
which is a file created by ``combine_hbonds.py``.  If more than one
file is given to the program, it will use only the first one.

Most of the options are the same as those for ``combine_hbonds.py``:
``-o`` (``--output-file``), ``-B`` (``--hbond-data-dir``), ``-y``
(``--occ-graph-only``), ``-R`` (``--resi-criteria``), ``-A``
(``--atom-criteria``), and ``-O`` (``--occ-thresh``).  Two options
differ: the sorting (``-s``, ``--sort``) and graph compression
(``-c``, ``--compress``) options.  The sorting option is
self-explanatory, but the other requires some explanation.  The graph
compression flag will shorten the length of the graph by combining
every pair of characters into a single character representing one
fifth (instead of one tenth) of the segment.  It should be noted that
the compressed graph is not as precise as the uncompressed graphs.
Each character in the compressed graph represents the occupancy
percentage over a pair of segments.  However, the exact occupancy
cannot always be determined solely from the ranges specified by the
characters from each segment.  For example, the actual occupancy of
the pair of segments represented by "x*" could correspond to either
"x" or "*", depending on the underlying percentages, which are
unknown.  In these cases, the character that is output is the one
corresponding to the most probable occupancy percentage.

compare_hbonds.py
-----------------

If H-bond data from multiple trajectories of the same system have been
processed by ``combine_hbonds.py``, the data from the different
trajectories can be compared with ``compare_hbonds.py``.

Documentation from the script::

  Usage: compare_hbonds.py FILE1 [ FILE2 [ ... ] ] [options] 
  FILE1 and additional optional files are files created by 
  combine_hbonds.py that contain datasets of H-bonds from different 
  trajectories of the same system.  The particular subset of the H-bonds 
  and the metric for sorting can be specified by the user.  
  
  Important notes:
  * Use the -i option to provide meaningful identifiers for the different
    trajectories.
  * The resi_criteria option takes a comma-separated string containing any
    or all of the following:
          - individual residue numbers
          - a range of numbers, separated by a '-'
          - strings associated with valid residue lists in the
              standard file 'residue_lists'
  * The atom_criteria option takes a comma-separated string containing
    the atom names to report.  In addition, the string can contain any
    of these strings:
          - 'bb_only': only H-bonds between two backbone atoms
          - 'not_bb': no H-bonds between two backbone atoms
          - 'protein_only': no H-bonds involving water
  * If the H-bonds are sorted by occ_pct (the occupancy percentage), any
    H-bond that has an occupancy greater than the occ_thresh value will
    be retained.  If the H-bonds are sorted by occ_diff (the difference
    between the largest and smallest occupancy percentages for the
    systems), donor, or acceptor, only those H-bonds with occ_diff
    greater than the occ_thresh value will be retained.
  * Note that ptraj uses a definition of H-bond donor and acceptor that is
    opposite of the normal convention.  This program follows the
    definitions of ptraj, in which the acceptor is the atom covalently 
    bonded to the hydrogen atom.
  
  Options:
    -h, --help            show this help message and exit
    -o OUTPUT_FILE, --output-file=OUTPUT_FILE
                          The name of the output file.  If None, the results
                          will be written to stdout.  [default: none]
    -B HBOND_DATA_DIR, --hbond-data-dir=HBOND_DATA_DIR
                          The directory that contains the H-bond data files.  If
                          None, the file names will be used without
                          modification, and output will be written to the
                          current directory.  [default: none]
    -i IDENTIFIERS, --identifiers=IDENTIFIERS
                          Comma-separated list of identifying strings for the
                          trajectories to be compared.  If None, the
                          trajectories will simply be numbered.  [default: none]
    -s SORT, --sort=SORT  The quantity used to sort the results.  Must be one of
                          "occ_diff" (occupancy difference), "occ_pct"
                          (occupancy percentage), "donor", or "acceptor".  The
                          occupancy difference is the difference between the
                          highest and lowest occupancy percentages for a
                          particular H-bond in the different trajectories.
                          [default: occ_diff]
    -c, --compress        Flag to compress the H-bond graph.  [default: False]
    -y, --occ-graph-only  Flag to report only the occupancy and graph data (no
                          distance or angle data).  [default: False]
    -R RESI_CRITERIA, --resi-criteria=RESI_CRITERIA
                          A comma- and dash-separated list of residue numbers to
                          include in the analysis.  [default: all]
    -A ATOM_CRITERIA, --atom-criteria=ATOM_CRITERIA
                          A comma-separated list of atom names to include in the
                          analysis.  [default: all]
    -O OCC_THRESH, --occ-thresh=OCC_THRESH
                          The minimum occupancy threshold that the H-bonds must
                          have to be reported.  [default: 0.0]


This program will allow the user to compare side-by-side the H-bond
data between two trajectories of the same system.  At least one input
file ``FILE1``, containing the H-bond data from ``combine_hbonds.py``,
is required.  In the output, the data for a given H-bond for all the
trajectories are shown next to each other.  The output has the same
format as shown in the example in the ``combine_hbond.py``
description, except that an identifying string is placed at the
beginning of each line to indicate the trajectory the data represents.
These identifiers can be input with the ``-i`` (``--identifiers``)
option.

::

  1w50A:  Glu165      N--H ... OD1  Asn162   36.00(  250) | ..xo.xxo-|-oo*o.....|
  1sgz :  Glu165      N--H ... OD1  Asn162   86.80(  250) |*o@@**x**@|*x@*x****@|


Most of the options are the same as those described for
``combine_hbonds.py`` and ``subset_hbonds.py``.  The only difference
is in the interaction of the sorting (``-s``, ``--sort``) and
occupancy threshold (``-O``, ``--occ-thresh``) options.  If the
H-bonds are sorted by occupancy difference, donor, or acceptor, the
output will include only those H-bonds for which the occupancy
*difference* is greater than the threshold.  If they are sorted by
occupancy percentage, the output will include only those H-bonds that
have at least one trajectory with an occupancy *percentage* greater
than the threshold.
