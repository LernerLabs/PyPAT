#!/usr/bin/env python

'''
Set up ptraj input files to run a series of H-bond calculations.
'''

import sys, os
from optparse import OptionParser
from pypat import hbond_tool_utils
parse_residue_list,parse_atom_list = hbond_tool_utils.parse_residue_list,hbond_tool_utils.parse_atom_list
from pypat import hbond_analysis_utils
compare_hbonds = hbond_analysis_utils.compare_hbonds

if __name__ == '__main__':

    usage = """%prog FILE1 [ FILE2 [ ... ] ] [options] 
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
  bonded to the hydrogen atom.""" 
    
    # Parse command line options

    parser = OptionParser(usage = usage)

    parser.add_option('-o', '--output-file',
           dest = 'output_file',
           help = 'The name of the output file.  If None, the results will be written to stdout.  [default: %default]')
    parser.add_option('-B', '--hbond-data-dir',
           dest = 'hbond_data_dir', 
	   help = 'The directory that contains the H-bond data files.  If None, the file names will be used without modification, and output will be written to the current directory.  [default: %default]')
    parser.add_option('-i', '--identifiers',
           dest = 'identifiers',
           help = 'Comma-separated list of identifying strings for the trajectories to be compared.  If None, the trajectories will simply be numbered.  [default: %default]')
    parser.add_option('-s', '--sort',
           dest = 'sort', default = 'occ_diff',
           choices = ['occ_diff', 'occ_pct', 'donor', 'acceptor'],
           help = 'The quantity used to sort the results.  Must be one of "occ_diff" (occupancy difference), "occ_pct" (occupancy percentage), "donor", or "acceptor".  The occupancy difference is the difference between the highest and lowest occupancy percentages for a particular H-bond in the different trajectories.  [default: %default]')
    parser.add_option('-c', '--compress', action = 'store_true',
           dest = 'compress', default = False,
           help = 'Flag to compress the H-bond graph.  [default: %default]')
    parser.add_option('-y', '--occ-graph-only', action = 'store_true',
           dest = 'occ_graph_only', default = False,
           help = 'Flag to report only the occupancy and graph data (no distance or angle data).  [default: %default]')
    parser.add_option('-R', '--resi-criteria',
           dest = 'resi_criteria', default = 'all',
           help = 'A comma- and dash-separated list of residue numbers to include in the analysis.  [default: %default]')
    parser.add_option('-A', '--atom-criteria',
           dest = 'atom_criteria', default = 'all',
           help = 'A comma-separated list of atom names to include in the analysis.  [default: %default]')
    parser.add_option('-O', '--occ-thresh', type = 'float',
           dest = 'occ_thresh', default = 0.0,
           help = 'The minimum occupancy threshold that the H-bonds must have to be reported.  [default: %default]')

    options, hbond_files = parser.parse_args()

    # Process options

    if options.identifiers:
	identifiers = options.identifiers.split(',')
    else:
	identifiers = []

    resi_criteria = parse_residue_list(options.resi_criteria, offset = 0)
    if not resi_criteria:
        sys.exit('ERROR:  No residues selected in residue string.\n')
    atom_criteria = parse_atom_list(options.atom_criteria)

    # Perform function

    compare_hbonds(hbond_files = hbond_files,
	identifiers = identifiers,
	output_file = options.output_file,
	resi_criteria = resi_criteria,
	atom_criteria = atom_criteria,
	occ_thresh = options.occ_thresh,
	occ_graph_only = options.occ_graph_only,
	sort = options.sort,
	compress = options.compress, 
	hbond_data_dir = options.hbond_data_dir  )

