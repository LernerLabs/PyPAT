#!/usr/bin/env python
import os,sys

import sys,os

if sys.platform == 'darwin':
    PYPAT_CODE_DIR = '/Users/mglerner/work/Dynamics-DHFR/'
elif sys.platform == 'linux2':
    PYPAT_CODE_DIR = '/users/mlerner/work/src/Dynamics-DHFR/'
    
sys.path.append(PYPAT_CODE_DIR)

if __name__ == '__main__':
    from pypat.runningptraj.write_ptraj_input_files import write_ptraj_input_files
    from pypat import tool_utils
    from optparse import OptionParser

    usage = tool_utils.usage + """
This script will emit the PDB file that you will use as a
reference structure later on.  It will live in <outputdir>/<structure>_ref.pdb.1
"""
    
    parser = OptionParser(option_class=tool_utils.MyOption,usage=usage)

    tool_utils.add_standard_options(parser)
    tool_utils.add_window_options(parser)

    parser.add_option("--input-dir",dest="inputdir",
                      default="./",
                      help="Directory that contains our input files.  It should contain the prmtop file and the mdcrd file.  [default: %default]")

    parser.add_option("--mdcrd",dest="mdcrd",
                      type="strlist",
                      default=["1sgz.nowat.0-6.mdcrd",],
                      help="Comma-separated list of mdcrd files")

    parser.add_option("--ps",dest="ps",
                      default=5,
                      type="int",
                      help="Number of ps per frame.  [default: %default]")

    parser.add_option("--align",dest="align",
                      default="all",
                      help="How to align. 'all' means 'rms first *'.  'none' means no alignment.  Any other string will be treated as the alignment string.  For example, if you say ':1-428@CA' the ptraj file will say 'rms first :1-428@CA'.  [default: %default]",
                      )
    parser.add_option("--strip-hydros",dest="strip_hydros",
                      default=False,
                      action="store_true",
                      help="Strip the hydrogens out during the ptraj runs.",
                      )
    parser.add_option("--strip-waters",dest="strip_waters",
                      default=False,
                      action="store_true",
                      help="Strip the waters out during the ptraj runs.  This assumes they're named WAT.",
                      )
    parser.add_option("--write-covar",dest="write_covar",
                      default=False,
                      action="store_true",
                      help="Write out the covariance matrix [default: %default]",
                      )
    parser.add_option("--other-ptraj-strips",dest="other_ptraj_strips",
                      default=[],
                      type="strlist",
                      help="Comma-separated list of other things that ptraj should strip.  For example, could be :WAT,:BOB and we would add two lines, one saying strip :WAT and one saying strip :BOB."
                      )
                      
    
    options,args = parser.parse_args()

    desired=tool_utils.get_desired(options)
    #
    # Build up the ptraj header.
    # 1) read in the mdcrd files
    # 2) strip things
    # 3) rms vs. first one
    #
    # VERY IMPORTANT NOTE: do the rms *after* you've done everything
    # else.  Otherwise, it tries to align things with the waters, etc.
    # and all of your correlations will be totally wrong.
    #
    ptraj_header = ''
    pdb_ptraj_header = ''
    pdb_ptraj_header += 'trajin %s 1 1 1\n\n'%os.path.join(options.inputdir,options.mdcrd[0])
    
    for mdcrd in options.mdcrd:
        ptraj_header += 'trajin %s\n'%os.path.join(options.inputdir,mdcrd)
    ptraj_header += '\n'

    if options.strip_hydros:
        ptraj_header += 'strip @H*\n'
        pdb_ptraj_header += 'strip @H*\n'
    if options.strip_waters:
        ptraj_header += 'strip :WAT\n'
        pdb_ptraj_header += 'strip :WAT\n'
    for strip in options.other_ptraj_strips:
        ptraj_header += 'strip %s\n'%strip
        pdb_ptraj_header += 'strip %s\n'%strip

    if options.align in 'none None NONE no NO No'.split():
        pass
    elif options.align == 'all':
        ptraj_header += 'rms first *\n'
        pdb_ptraj_header += 'rms first *\n'
    else:
        ptraj_header += 'rms first %s\n'%options.align
        pdb_ptraj_header += 'rms first %s\n'%options.align
            
    
        
        
    write_ptraj_input_files(dir_containing_ptraj_files=options.outputdir,
                            ptraj_output_dir=os.path.join(options.outputdir,options.structurename),
                            filename_prefix=options.structurename,
                            desired=desired,
                            ps_per_frame=options.ps,
                            ptraj_header=ptraj_header,
                            fname_template='calculate_'+options.structurename+'_%s_correl_and_covar.ptraj',
                            write_out_pdb_file=os.path.join(options.outputdir,options.structurename+'_ref.pdb'),
                            write_covar=options.write_covar,
                            pdb_ptraj_header=pdb_ptraj_header,
                            )
    
    
