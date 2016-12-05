Bridging-water analysis
=======================

Executive Summary
-----------------

An understanding of bridging-water molecules is critical in the study
of structure, dynamics, and function of proteins and nucleic
acids. While several programs (such as ptraj) allow for ane analysis
of hydrogen bonds, the analysis of bridging waters is significantly
more tedius. Typically, users must extract this information either
through detailed examination of hydrogen bonds between pairs of
protein atoms and water atoms, or through visual examination of MD
trajectories.

PyPAT greatly simplifies this process, providing two scripts
(``collect_water_bridges.py`` and
``display_bridging_interactions.py``) that extract and analyze
bridging interactions throughout an MD trajectory.

A typical invocation for a trajectory and topology named ``mysim.trj``
and ``mysim.top`` might look like::

  prompt$ pymol -qcr collect_water_bridges.py -- --name=mysim
  prompt$ display_bridging_interactions.py --name=mysim --resi-criteria 5,6,10-24

These scripts are able to analyze proteins, nucleic acids, and
custom-defined ligands. In order to simplify the wording below, all of
these will be referred to as "protein".

Phase 1 (collect_water_bridges.py)
----------------------------------

Description
~~~~~~~~~~~

The first step is the most time consuming: extracting the bridging
information from the MD simulation. Using PyMOL as a backend, all
water molecules for which the Oxygen is within 4.0 Angstroms of the
protein are examined. Information about the distance and angle of each
interaction between these molecules and the protein is recorded. This
information can be further refined in the next step, so we recommend
using loose constraints; this helps avoid the need to rerun
``collect_water_bridges.py``. 

We note again that this script must be run through PyMOL::

  pymol -qcr path/to/collect_water_bridges.py -- --name=mysimulation

etc. The ``--`` after the script name is required. 

The default options for distance and angle cutoffs are usually
correct. Users must explicitly specifiy the number of steps in the MD
trajectory (``--num-steps``), as well as the name of the trajectory
(``--name``). The trajectory and topology file must be named
``<name>.top`` and ``<name>.trj`` respectively.

PyMOL slows down if it processes too many frames at once. Therefore,
the trajectory is analyzed in chunks of 500 frames at a time. We find
this to be a generally useful chunk-size, but users can change it via
the ``chunk-size`` command-line option.

Defining new ligands
~~~~~~~~~~~~~~~~~~~~

The script comes with definitions for standard protein and nucleic
acid hydrogen-bonding interactions. Users may wish to change these, or
to add definitions for new ligands. This is easily accomplished by
editing the file ``hbond_definitions.py`` (installed under
``pypat/hbond/`` when the scripts are installed). As an example, here
is the function that defines donors and acceptors for NADPH::

  def select_nap_donors_and_acceptors():
      """
      Ligand specific selections for NADPH (NAP)
      """
      #-- NADPH
      #acceptor mask :NAP@N6A  :NAP@H61
      cmd.select("prot_donors","prot_donors or (resn NAP and name H61)")
      #acceptor mask :NAP@N6A  :NAP@H62
      cmd.select("prot_donors","prot_donors or (resn NAP and name H62)")
      #acceptor mask :NAP@O'A3 :NAP@HOA3
      cmd.select("prot_donors","prot_donors or (resn NAP and name HOA3)")
      #acceptor mask :NAP@O'N3 :NAP@HON3
      cmd.select("prot_donors","prot_donors or (resn NAP and name HON3)")
      #acceptor mask :NAP@O'N2 :NAP@HON2
      cmd.select("prot_donors","prot_donors or (resn NAP and name HON2)")
      #acceptor mask :NAP@N7N  :NAP@H72
      cmd.select("prot_donors","prot_donors or (resn NAP and name H72)")
      #acceptor mask :NAP@N7N  :NAP@H71
      cmd.select("prot_donors","prot_donors or (resn NAP and name H71)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name N1A)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name N3A)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name N7A)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OA23)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OA22)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OA24)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'A2)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'A3)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'A4)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'A5)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OPA1)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OPA2)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OPN1)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O3P)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OPN2)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'N5)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'N4)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'N3)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'N2)")
      cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O7N)")

After defining such a function, it must be added to the
``do_standard_selections`` function at the bottom of ``hbond_definitions.py``.


Documentation from the script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  Usage: 
      Run this like:
  
      pymol -qcr collect_water_bridges.py -- --name=myprotein --dist-cutoff=3.5
  
      Do not forget the double dashes after the script name.
      
  
  Options:
    -h, --help            show this help message and exit
    -n NAME, --name=NAME  Trajectory and topology must be named name.trj and
                          name.top respectively. [default: nrna]
    -c CHUNKSIZE, --chunk-size=CHUNKSIZE
                          How many MD steps to process at a time. If you do too
                          many at a time, PyMOL will slow down. Too few, and
                          you're wasting time starting/stopping PyMOL. [default:
                          500]
    -s NUMSTEPS, --num-steps=NUMSTEPS
                          number of steps in your MD trajectory. [default: 5000]
    -d DISTCUTOFF, --dist-cutoff=DISTCUTOFF
                          Heavy atom to heavy atom distance cutoff. [default:
                          4.0]
    -a ANGLECUTOFF, --angle-cutoff=ANGLECUTOFF
                          Angle cutoff. If heavy:hydro:heavy angle must be
                          greater than this. [default: 0.0]

Phase 2 (display_bridging_interactions.py)
------------------------------------------

Description
~~~~~~~~~~~

Phase 1 records data in a water-centric format. That is, interactions
are recorded and described for each water molecule. Phase 2 inverts
this, and displays bridging interactions as protein-water-protein
triplets. There are several subtleties involved in this process, and
users are strongly advised to read the paper. Among the relevant
options are:

minrequireddwelltime 
  Bridging interactions that do not persist for at
  least this long are ignored

looseness
  This allows for gaps in the occupancy of a bridging interaction.
  Suppose a bridging interaction is present for 300 picoseconds,
  absent for 2 picoseconds, and present for another 200 picoseconds.
  If the looseness is greater than or equal to 2 picoseconds, this 
  will be treated as a single 502 picosecond interaction.

minocc
  Bridging interactions that are absent for substantial percentages of the
  trajectory can be automatically filtered out.

resi-criteria
  A complete list of bridging interactions will contain an overwhelming 
  amount of information, including many surface interactions in regions
  that may not be interesting to the user. With this option, a user can
  specify a list of residues of interest. Bridging interactions that do not
  involve at least one of these residues are ignored. The input format is
  fairly general (e.g. "5,6,7-12"). Users who are comfortable with Python
  can find examples in the code (``tool_utils.py``) of how to define residue
  groups such as "loop A", etc.


Documentation from the script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  Usage: display_bridging_interactions.py [options]
  
  Options:
    -h, --help            show this help message and exit
    -n NAME, --name=NAME  Trajectory and topology must be named name.trj and
                          name.top respectively. [default: nrna]
    -t TIMESTEP, --timestep=TIMESTEP
                          trajectory timestep in picoseconds. [default: 5]
    -m MINREQUIREDDWELLTIME, --min-dwell-time=MINREQUIREDDWELLTIME
                          minimum required dwell time. [default: 3]
    -l LOOSENESS, --looseness=LOOSENESS
                          looseness. [default: 2]
    -d DISTCUTOFF, --dist-cutoff=DISTCUTOFF
                          Heavy atom to heavy atom distance cutoff. [default:
                          3.5]
    -a ANGLECUTOFF, --angle-cutoff=ANGLECUTOFF
                          Angle cutoff. If heavy:hydro:heavy angle must be
                          greater than this. [default: 0.0]
    -o MINOCC, --min-occ=MINOCC
                          Minimum percentage of the trajectory for which this
                          interaction must be occupied. [default: 0.4]
    -r DIR, --dir=DIR     Directory in which the name_hbond_*_*.txt files are
                          located. [default: .]
    -R RESI_CRITERIA, --resi-criteria=RESI_CRITERIA
                          Restrict the output to BWIs where at least one side
                          involves a residue in this list. [default: none]

