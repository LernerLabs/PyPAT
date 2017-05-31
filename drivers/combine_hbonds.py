#!/cluster/home2/mglerner/anaconda3/bin/python
#!/usr/bin/env python

'''
Combine hbond output from a series of ptraj runs into one data set.
'''

import sys, os
from optparse import OptionParser
from pypat import md_analysis_utils
get_resinum_to_resi_map = md_analysis_utils.get_resinum_to_resi_map
from pypat import hbond_tool_utils
parse_residue_list,parse_atom_list = hbond_tool_utils.parse_residue_list,hbond_tool_utils.parse_atom_list
from pypat import hbond_analysis_utils
combine_hbonds = hbond_analysis_utils.combine_hbonds

if __name__ == '__main__':

    usage = """%prog FILE1 [ FILE2 [ ... ] ] [options] 
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
	- 'protein_only': no H-bonds involving water""" 
    
    # Parse command line options

    parser = OptionParser(usage = usage)

    parser.add_option('-o', '--output-file',
           dest = 'output_file', 
           help = 'The name of the output file.  If None, the results will be written to stdout.  Supplying a name is recommended.  [default: %default]')
    parser.add_option('-B', '--hbond-data-dir',
           dest = 'hbond_data_dir', 
           help = 'The directory that contains the H-bond data files.  If None, the file names will be used without modification, and output will be written to the current directory.  [default: %default]')
    parser.add_option('-g', '--segment-size', type = 'int',
           dest = 'segment_size', default = 1000,
           help = 'The number of frames included in each segment of the trajectory.  [default: %default]')
    parser.add_option('-p', '--prmtop-file',
           dest = 'prmtop_file', 
           help = 'Amber parameter/topology file or PDB file for the system.  [default: %default]')
    parser.add_option('-D', '--prmtop-dir',
           dest = 'prmtop_dir', default = None,
           help = 'The directory that contains the prmtop (or PDB) file.  If None, the name of the file will be used without modification.  [default: %default]')
    parser.add_option('-r', '--resi-offset', type = 'int',
           dest = 'resi_offset', default = 0,
           help = 'The offset between the residue numbers in the prmtop (or PDB) file and the actual residue numbers.  [default: %default]')
    parser.add_option('-a', '--aa-code', type = 'int',
	   dest = 'aa_code', default = 3,
	   help = 'Must be 1 or 3.  Indicates the use of 1- or 3-letter amino acid codes in residue names.  [default: %default]')
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

    if options.aa_code != 1 and options.aa_code != 3:
	print('Warning:  Illegal amino acid code.  Must be 1 or 3.\n' + \
	      '  Will use the default of 3.')
	options.aa_code = 3

    if options.prmtop_dir != None:
	options.prmtop_file = os.path.join(options.prmtop_dir, options.prmtop_file)
    resi_map = get_resinum_to_resi_map(options.prmtop_file, 
		offset = options.resi_offset, aa_code = options.aa_code)

    resi_criteria = parse_residue_list(options.resi_criteria, offset = 0)
    if not resi_criteria:
	sys.exit('ERROR:  No residues selected in residue string.\n')
    atom_criteria = parse_atom_list(options.atom_criteria)

    # Perform function

    combine_hbonds(hbond_files = hbond_files,
	segment_size = options.segment_size,
	resi_map = resi_map,
	output_file = options.output_file,
	resi_criteria = resi_criteria,
	atom_criteria = atom_criteria,
	occ_thresh = options.occ_thresh,
	occ_graph_only = options.occ_graph_only,
	hbond_data_dir = options.hbond_data_dir  )

