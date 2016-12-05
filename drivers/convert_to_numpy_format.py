#!/usr/bin/env python

import glob,bz2
import scipy
import scipy.io
import os


if __name__ == '__main__':
    from optparse import OptionParser
    usage = """This will convert the .dat or .dat.bz2 files to numpy versions.
    It will not automatically delete the .dat(.bz2) files.  If your input
    files are bz2, your output files will be too.

    If you already have a corresponding .numpy or .numpy.bz2 file, we won't
    write out a new file.
    """
    parser = OptionParser(usage=usage)

    parser.add_option('--dir',dest='dir',default='.',
                      help="Directory in which the files reside. [default: %default]")
    parser.add_option("--structure-name",dest="structurename",
                      default="1rx1",
                      help="Name of your structure.  E.g. 1RX1 or 1SGZ")
    parser.add_option('--compression',dest='compression',
                      default='',
                      help="Type of compression currently used on files.  Leave blank for uncompressed, .bz2 if they're .bz2 files. Note that it's '.bz2' not 'bz2'. [default: %default]")
    parser.add_option('--all-dat-files',dest='justbig',
                      default=True,
                      action='store_false',
                      help="By default, we will only convert the all_atom_correlmat files.  If you use this option, we will convert all dat files.  Don't forget, though, that you'll still have to call this command twice if you have some files that are bzipped and some that are not.",
                      )

    options,args = parser.parse_args()
    if options.justbig:
        pattern = os.path.join(options.dir,options.structurename+'_NS*_all_atom_correlmat*.dat'+options.compression)
    else:
        pattern = os.path.join(options.dir,options.structurename+'_NS*.dat'+options.compression)
    fnames = glob.glob(pattern)
    print("I will convert these files",pattern,":",fnames)
    for fname in fnames:
        if fname.endswith('bz2'):
            head = os.path.splitext(fname)[0]
        else:
            head = fname
        numpy_ver = head+'.numpy'
        numpybz_ver = head+'.numpy.bz2'
        if os.path.exists(numpy_ver) or os.path.exists(numpybz_ver):
            print("Skipping",fname)
            continue
        if fname.endswith('bz2'):
            print("Doing",fname,numpybz_ver)
            data=scipy.io.read_array(bz2.BZ2File(fname))
            data.tofile(numpy_ver)
            os.system('bzip2 %s'%numpy_ver)
        else:
            print("Doing",fname,numpy_ver)
            data=scipy.io.read_array(file(fname))
            data.tofile(numpy_ver)
        del data

    
