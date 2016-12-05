#!/usr/bin/env python

'''
Set up ptraj input files to run a series of H-bond calculations.
'''

import sys, os
from optparse import OptionParser

if __name__ == '__main__':

    usage = """%prog [options]
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
  documentation."""
    
    # Parse command line options

    parser = OptionParser(usage = usage)

    parser.add_option('-x', '--mdcrd-file',
           dest = 'mdcrd_file',
           help = 'Coordinate file to use for H-bond analysis (REQUIRED).  [default: %default]')
    parser.add_option('-p', '--prmtop-file',
           dest = 'prmtop_file', 
           help = 'Amber parameter/topology file to use for H-bond analysis.  If None, the name will be guessed by replacing the extension of the coordinate file name with ".prmtop"  [default: %default]')
    parser.add_option('-D', '--mdcrd-prmtop-dir',
           dest = 'mdcrd_prmtop_dir', 
           help = 'The directory that contains the coordinate and prmtop file.  If None, the names of the files will be used without modification.  [default: %default]')
    parser.add_option('-o', '--output-file-base',
           dest = 'output_file_base', default = 'segment',
           help = 'The first part of the name of the ptraj input files that will be created.  The files will be named OUTPUT_FILE_BASEnn.in, where nn is the two-digit segment number.  Once ptraj is run (with run_hbond_ptraj), the output files containing the H-bond data will be named OUTPUT_FILE_BASEnn.out.  [default: %default]')
    parser.add_option('-B', '--hbond-data-dir',
           dest = 'hbond_data_dir', default = '.',
           help = 'The directory where the ptraj input files will be placed.  [default: %default]')
    parser.add_option('-b', '--begin-frame', type = 'int',
           dest = 'begin_frame', default = 1,
           help = 'The number of the first frame to use in the H-bond analysis.  [default: %default]')
    parser.add_option('-e', '--end-frame', type = 'int',
           dest = 'end_frame', default = 10000,
           help = 'The number of the last frame to use in the H-bond analysis.  [default: %default]')
    parser.add_option('-g', '--segment-size', type = 'int',
           dest = 'segment_size', default = 1000,
           help = 'The number of frames included in each segment of the trajectory.  [default: %default]')
    parser.add_option('-n', '--num-resi', type = 'int',
           dest = 'num_resi', 
           help = 'The number of residues in the protein.  [default: 10000]')
    parser.add_option('-d', '--dist-cutoff', type = 'float',
           dest = 'dist_cutoff', default = 3.0,
           help = 'The distance cutoff that defines whether or not an interaction is an H-bond.  [default: %default]')
    parser.add_option('-a', '--angle-cutoff', type = 'float',
           dest = 'angle_cutoff', default = 120.0,
           help = 'The angle cutoff that defines whether or not an interaction is an H-bond.  [default: %default]')
    parser.add_option('-s', '--no-self', action = 'store_false',
           dest = 'self', default = True,
           help = 'Flag to turn off the inclusion of H-bonds between atoms within the same residue.  [default: %default]')
    parser.add_option('-S', '--solvent', action = 'store_true',
           dest = 'solvent', default = False,
           help = 'Flag to include solvent-protein interactions in the analysis.  [default: %default]')
    parser.add_option('-m', '--mask-file',
           dest = 'mask_file', default = None,
           help = 'File that contains extra lines to include in the ptraj input files, primarily to include masks for ligands.  [default: mask]')

    options, args = parser.parse_args()

    mdcrd_file = options.mdcrd_file
    prmtop_file = options.prmtop_file
    mdcrd_prmtop_dir = options.mdcrd_prmtop_dir
    output_file_base = options.output_file_base
    begin_frame = options.begin_frame
    end_frame = options.end_frame
    segment_size = options.segment_size
    num_resi = options.num_resi
    dist_cutoff = options.dist_cutoff
    angle_cutoff = options.angle_cutoff
    self = options.self
    solvent = options.solvent
    mask_file = options.mask_file
    hbond_data_dir = options.hbond_data_dir

    # Do error check for file names

    if not mdcrd_file:
	sys.exit('ERROR:  No coordinate file given (-x, --mdcrd-file option).  For help, \n' + \
                 '  use the -h option.\n')
    if not prmtop_file:
        prmtop_file = os.path.splitext(options.mdcrd_file)[0] + '.prmtop'
        print 'Warning:  No prmtop file specified.  Will use ' + prmtop_file + '.' 
    if mdcrd_prmtop_dir:
        mdcrd_file = os.path.join(mdcrd_prmtop_dir, mdcrd_file)
        prmtop_file = os.path.join(mdcrd_prmtop_dir, prmtop_file)
    
    if not os.path.exists(mdcrd_file):
        sys.exit('ERROR:  Coordinate file ' + mdcrd_file + '\n' + \
                 '  does not exist.')
    if not os.path.exists(prmtop_file):
        sys.exit('ERROR:  prmtop file ' + prmtop_file + '\n' + \
                 '  does not exist.')

    # Read through prmtop file to determine which residue types
    # are present and also get the residue numbers of the prolines

    prmtop_f = file(prmtop_file)

    prolines = []    # list of proline residue numbers
    resi_num = 1
    is_present = {}  # dict of residue types present in the file
    residue_section = False
    for line in prmtop_f:
	if line.startswith('%FLAG RESIDUE_POINTER'):
	    break
	if line.startswith('%FLAG RESIDUE_LABEL'):
	    residue_section = True
	if not residue_section or line.startswith('%F'):
	    continue
	else:
	    residue_names = line.split()
	    for resi_name in residue_names:
		is_present[resi_name] = True
		if resi_name == 'PRO' and resi_num > 1:
		    prolines.append(resi_num)
		resi_num += 1

    resi_types = '''ALA ARG ASH ASN ASP CYM CYS CYX GLH GLN GLU GLY
		    HID HIE HIP ILE LEU LYN LYS MET PHE PRO SER THR
		    TRP TYR VAL'''.split()

    for resi_type in is_present:
	if resi_type not in resi_types:
		if not (solvent and resi_type == 'WAT'):
		    print 'Warning:  Unrecognized residue type ' + resi_type + ' is present in system.  Atoms in\n' + \
			  '  this residue will not be explicitly included in H-bond analysis unless its\n' + \
			  '  donors and acceptors are specified in the mask file.'

    # Set non-None defaults for the options that need them

    if mask_file == None:
	mask_file = 'mask'
	mask_warning = ''
    else:
	mask_warning = 'Warning:  Could not open ' + mask_file + '.  Ignoring -m option.'
    extra_masks = ''	
    try:
	mask_f = file(mask_file)
    except IOError:
	if mask_warning:
	    print mask_warning
    else:
	for line in mask_f:
	    extra_masks += line

    if num_resi == None:
	print 'Warning:  No num_resi option (-n, --num-resi) option specified.  Will use\n' + \
	      '  default value of 10000.  If your protein has fewer than 10000 residues\n' + \
              '  and there is solvent in the system, some of the solvent oxygen atoms\n' + \
              '  will be treated explicitly as donors.'
	num_resi = 10000

    # Create string of nonproline residues (- and , separated
    # like in normal ptraj input) based on proline list

    for resinum in range(2, num_resi + 1):
	if resinum in prolines:
            prolines = prolines[1:]
        else:
            prolines = [resinum - 1] + prolines + [num_resi + 1]
            break
    nonpro_str = ''
    for i in range(len(prolines) - 1):
	diff = prolines[i + 1] - prolines[i]
	if diff > 1:
            nonpro_str += str(prolines[i] + 1)
	    if diff > 2:
		nonpro_str += '-' + str(prolines[i + 1] - 1)
	    if i != len(prolines) - 2:
		nonpro_str += ','

    # Set strings of H-bond donor and acceptor atoms

    donor_mask = {}
    donor_mask['ASH'] = 'donor mask :ASH@OD1\n' + \
			'donor mask :ASH@OD2'
    donor_mask['ASN'] = 'donor mask :ASN@OD1'
    donor_mask['ASP'] = 'donor mask :ASP@OD1\n' + \
			'donor mask :ASP@OD2'
    donor_mask['CYM'] = 'donor mask :CYM@SG'
    donor_mask['GLH'] = 'donor mask :GLH@OE1\n' + \
			'donor mask :GLH@OE2'
    donor_mask['GLN'] = 'donor mask :GLN@OE1'
    donor_mask['GLU'] = 'donor mask :GLU@OE1\n' + \
			'donor mask :GLU@OE2'
    donor_mask['HID'] = 'donor mask :HID@NE2'
    donor_mask['HIE'] = 'donor mask :HIE@ND1'
    donor_mask['LYN'] = 'donor mask :LYN@NZ'
    donor_mask['SER'] = 'donor mask :SER@OG'
    donor_mask['THR'] = 'donor mask :THR@OG1'
    donor_mask['TYR'] = 'donor mask :TYR@OH'

    acceptor_mask = {}
    acceptor_mask['ARG'] = 'acceptor mask  :ARG@NH2 :ARG@HH21\n' + \
			   'acceptor mask  :ARG@NH2 :ARG@HH22\n' + \
			   'acceptor mask  :ARG@NH1 :ARG@HH11\n' + \
			   'acceptor mask  :ARG@NH1 :ARG@HH12\n' + \
			   'acceptor mask  :ARG@NE  :ARG@HE'
    acceptor_mask['ASH'] = 'acceptor mask  :ASH@OD2 :ASH@HD2'
    acceptor_mask['ASN'] = 'acceptor mask  :ASN@ND2 :ASN@HD21\n' + \
			   'acceptor mask  :ASN@ND2 :ASN@HD22'
    acceptor_mask['GLH'] = 'acceptor mask  :GLH@OE2 :GLH@HE2'
    acceptor_mask['GLN'] = 'acceptor mask  :GLN@NE2 :GLN@HE21\n' + \
			   'acceptor mask  :GLN@NE2 :GLN@HE22'
    acceptor_mask['HID'] = 'acceptor mask  :HID@ND1 :HID@HD1'
    acceptor_mask['HIE'] = 'acceptor mask  :HIE@NE2 :HIE@HE2'
    acceptor_mask['HIP'] = 'acceptor mask  :HIP@ND1 :HIP@HD1\n' + \
			   'acceptor mask  :HIP@NE2 :HIP@HE2'
    acceptor_mask['LYN'] = 'acceptor mask  :LYN@NZ  :LYN@HZ2\n' + \
			   'acceptor mask  :LYN@NZ  :LYN@HZ3'
    acceptor_mask['LYS'] = 'acceptor mask  :LYS@NZ  :LYS@HZ1\n' + \
			   'acceptor mask  :LYS@NZ  :LYS@HZ2\n' + \
			   'acceptor mask  :LYS@NZ  :LYS@HZ3'
    acceptor_mask['SER'] = 'acceptor mask  :SER@OG  :SER@HG'
    acceptor_mask['THR'] = 'acceptor mask  :THR@OG1 :THR@HG1'
    acceptor_mask['TRP'] = 'acceptor mask  :TRP@NE1 :TRP@HE1'
    acceptor_mask['TYR'] = 'acceptor mask  :TYR@OH  :TYR@HH'

    # Create the ptraj input files and run_all file

    try:
        run_file = os.path.join(hbond_data_dir, 'run_hbond_ptraj')
        run_f = file(run_file, 'w')
    except IOError:
        sys.exit('ERROR:  Could not open ' + run_file + '.\n')

    num_segments = (end_frame - begin_frame + 1) / segment_size
    remainder = (end_frame - begin_frame + 1) % segment_size
    if remainder:
	print 'Warning:  The number of frames is not a multiple of the segment size.\n' + \
              '  Ignoring final ' + str(remainder) + ' frames.'
        end_frame -= remainder

    for segment_index in range(num_segments):

        # Write ptraj input file

	filename_base = output_file_base
	if segment_index < 9:
	    filename_base += '0'
	filename_base += str(segment_index + 1)
        output = ''

	# Trajectory to read in

        output += 'trajin ' + mdcrd_file + ' ' + \
                   str(begin_frame + segment_index*segment_size) + ' ' + \
                   str(begin_frame + (segment_index+1)*segment_size - 1) + '\n'

	# List the donor atoms of sidechains

	output += '\n# List of potential H-bond donors\n'
	for resi_type in is_present:
	    if donor_mask.has_key(resi_type):
		output += donor_mask[resi_type] + '\n'

	# List the acceptor atoms of sidechains

	output += '\n# List of potential H-bond acceptors\n'
	for resi_type in is_present:
	    if acceptor_mask.has_key(resi_type):
		output += acceptor_mask[resi_type] + '\n'

	# List the backbone atoms

	output += '\n#-- Backbone donors and acceptors for ' + \
				'this particular molecule\n' + \
		    '#   N-H for prolines do not exist so ' + \
				'are not in the mask\n'
        output += 'donor mask :1-' + str(num_resi) + '@O\n' + \
                    'acceptor mask  :' + \
                    nonpro_str + \
                    '@N :2-' + str(num_resi) + '@H'

	output += '\n#-- Terminal residues have different atom names\n' + \
		    'donor mask @OXT\n' + \
		    'acceptor mask :1@N :1@H1\n' + \
		    'acceptor mask :1@N :1@H2\n' + \
		    'acceptor mask :1@N :1@H3\n'

	# List the extra lines from mask file

	output += '\n#-- Masks supplied by user\n' + extra_masks
	if not extra_masks:
	    output += '#   None given\n'

	# hbond process line

        output += '\nhbond distance ' + str(dist_cutoff) + \
                    ' angle ' + str(angle_cutoff)
        if solvent:
            output += ' solventneighbor 6 solventdonor WAT O' + \
		      ' solventacceptor WAT O H1 solventacceptor WAT O H2'
        output += ' series hb_all out ' + filename_base + '.out'
        if self:
            output += ' includeself\n'
        else:
            output += '\n'

	# Output the ptraj file for the segment

        try:
            output_file = filename_base + '.in'
	    output_file = os.path.join(hbond_data_dir, output_file)
            output_f = file(output_file, 'w')
        except IOError:
            print 'ERROR:  Could not open ' + output_file + '.\n'
            sys.exit('  Writing out ptraj input file to stdout:\n' + output)
        else:
            output_f.write(output + '\n')
            output_f.close()

        # Update runfile

        run_file_output = 'ptraj ' + prmtop_file + \
			  ' ' + os.path.split(output_file)[-1] + '\n'
        run_f.write(run_file_output)

    run_f.close()
    os.system('chmod +x ' + run_file)

