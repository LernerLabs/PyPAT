#!/cluster/home2/mglerner/anaconda3/bin/python
#!/usr/bin/env python

import sys,os

if __name__ == '__main__':
    from pypat.md_analysis_utils import write_max_min_ca_resi_versions
    from pypat import tool_utils
    from optparse import OptionParser    
    #
    # Standard procedure for writing out the specific correlation and covariance
    # matrices that we want.
    #

    parser = OptionParser(option_class=tool_utils.MyOption,usage=tool_utils.usage)


    tool_utils.add_standard_options(parser)
    tool_utils.add_window_options(parser)
    
    parser.add_option("--non-ca-resis",dest="non_ca_resis",type="zerobasedintlist",
                      default=[],
                      help="Comma separated list of residues that don't contain alpha carbons.  We need this to make some of our output images.  [default: %default], but you could say 160,161 for example")


    options,args = parser.parse_args()

    if options.singletime is None:
        desired=tool_utils.get_desired(options)
        for times in [i[-1] for i in desired]:
            sys.stdout.write('doing %s %s'%(times,options.structurename))
            sys.stdout.flush()
            write_max_min_ca_resi_versions(os.path.join(options.outputdir,options.structurename,'%s_%s_all_atom_correlmat.dat'%(options.structurename,times)),
                                           os.path.join(options.outputdir,options.structurename+'_ref.pdb.1'),
                                           non_ca_resis=[i+1 for i in options.non_ca_resis],
                                           )
    else:
            write_max_min_ca_resi_versions(os.path.join(options.outputdir,options.structurename,'%s_%s_all_atom_correlmat.dat'%(options.structurename,options.singletime)),
                                           os.path.join(options.outputdir,options.structurename+'_ref.pdb.1'),
                                           non_ca_resis=[i+1 for i in options.non_ca_resis],
                                           )
