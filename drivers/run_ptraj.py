#!/usr/bin/env python

import sys,os

if sys.platform == 'darwin':
    PYPAT_CODE_DIR = '/Users/mglerner/work/Dynamics-DHFR/'
elif sys.platform == 'linux2':
    PYPAT_CODE_DIR = '/users/mlerner/work/src/Dynamics-DHFR/'
    
sys.path.append(PYPAT_CODE_DIR)

if __name__ == '__main__':
    from pypat import tool_utils
    import glob
    from optparse import OptionParser

    usage = tool_utils.usage + """
Please also make sure that ptraj is installed and in your path.

This mimics the following shell script:

    export PTRAJ_INPUT_DIR=/users/mlerner/tmp/spronk/ptraj_files/
    export PRMTOP=/users/spronk/temp/1sgz/1sgz.nowat.prmtop
    for f in $PTRAJ_INPUT_DIR/*.ptraj; do ptraj $PRMTOP <$f > $f.out; done
"""
    parser = OptionParser(usage=tool_utils.usage)

    tool_utils.add_standard_options(parser)

    parser.add_option("--input-dir",'-i',dest="inputdir",
                      default="./",
                      help="Directory that contains our input files.  It should contain the prmtop file and the mdcrd file.  [default: %default]")
    parser.add_option("--prmtop",'-p',dest="prmtop",
                      default="1sgz.nowat.prmtop",
                      help="Name of prmtop file")

    options,args = parser.parse_args()

    filenames = glob.glob(os.path.join(options.outputdir,'calculate_%s*_correl_and_covar.ptraj'%options.structurename))
    filenames += glob.glob(os.path.join(options.outputdir,'write_%s_ref_pdb.ptraj'%options.structurename))
    print "I found these ptraj files to run:",filenames
    for f in filenames:
        prog,args = 'ptraj',(os.path.join(options.inputdir,options.prmtop), f)
        retcode,progout = tool_utils.run(prog,args)
        outf = file(f+'.out','w')
        outf.write(progout)
        outf.close()
        print "Ran",prog,args,"with result",retcode

