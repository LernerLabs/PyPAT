#!/usr/bin/env python

"""

Various utility code that gets used throughout the command-line tools.

Recent changes

 - Moved scipy and numpy imports into relevant functions so that
   other python builds (e.g. PyMOL) can easily use this.

"""

from __future__ import division
import os,sys,bz2
from copy import copy
from optparse import Option, OptionValueError

usage = """usage: %prog [options]

Please make sure that you have created the following directories:

  output-dir
  output-dir/structure-name
  output-dir/images
"""

############################
#
# New Optparse option types
#
############################

# Define a new option type, a comma-separated list of integers.
def check_zerobasedintlist(option,opt,value):
    """
    This automatically converts one-based indices to zero-based
    indices.  If you don't want that, make sure to convert back
    yourself!!

    It also knows about all of our lists of residues.
    """
    #
    # NOTE: if you want residues 9-24, give
    #       range(8,24) because we're zero-indexed.
    #       
    standard_residue_lists = {
        # Loops
        'dhfr_fg':range(115,132),
        'dhfr_cd':range(63,71),
        'dhfr_m20':range(8,24),
        'dhfr_gh':range(141,150),
        'dhfr_subdomain1':range(37) + range(106,159),
        'dhfr_subdomain2':range(38,106),
        # Beta strands
        'dhfr_BA':[i-1 for i in range(  2,  5+1)],
        'dhfr_BB':[i-1 for i in range( 39, 43+1)], #
        'dhfr_BC':[i-1 for i in range( 58, 62+1)], #
        'dhfr_BD':[i-1 for i in range( 73, 75+1)], #
        'dhfr_BE':[i-1 for i in range( 91, 95+1)],
        'dhfr_BF':[i-1 for i in range(107,115+1)],
        'dhfr_BG':[i-1 for i in range(133,135+1)],
        'dhfr_BH':[i-1 for i in range(151,158+1)],
        # Alpha helices
        'dhfr_AB':[i-1 for i in range( 25, 35+1)],
        'dhfr_AC':[i-1 for i in range( 44, 50+1)], #
        'dhfr_AE':[i-1 for i in range( 78, 85+1)], #
        'dhfr_AF':[i-1 for i in range( 97,106+1)],

        # Networks
        'dhfr_AgarwalJPhysChemBNetwork':[i-1 for i in (7,14,15,27,28,31,40,41,43,44,46,47,54,61,62,63,100,113,122)],
        'dhfr_WatneyHammes-SchifferJPhysChemBNetwork_shared':[i-1 for i in (28,42,44,50,51,52,53,64,77,90)],
        'dhfr_WatneyHammes-SchifferJPhysChemBNetwork_ecolionly':[i-1 for i in (19,20,45,47,61,62,63,67,68,72,74,98,99,102,108,129,149,155)],
        # Mutants
        'dhfr_rb_catcorr_mutants':[i-1 for i in (9,44,45,46,54,121,122)],
        'dhfr_rb_catnoncorr_mutants':[i-1 for i in (28,31,100,113)],
        'dhfr_rb_noncatnoncorr_mutants':[i-1 for i in (85,88,137,145,152,153,155)],
        'dhfr_rb_noncatcorr_mutants':[i-1 for i in (49,67,14,22)],
        # Other DHFR
        'dhfr_allosteric_site':[i-1 for i in        (      26,   29,30,   33,   111,        137,139,141,153,    155)],
        'dhfr_moe_allosteric_exposed':[i-1 for i in (      26,   29,30,   33,               137,139,141,153,    155)],
        'dhfr_moe_allosteric_tunnel': [i-1 for i in (5,6,7,   27,   30,31,   34,111,112,113,            153,154)],
        'dhfr_allosteric_everything': [i-1 for i in (5,6,7,26,27,29,30,31,33,34,111,112,113,137,139,141,153,154,155)],

        'dhfr_allosteric_nonexposed': [i-1 for i in (5,6,7,   27,      31,   34,    112,113,                154)],

        'dhfr_mtx_contact':[i-1 for i in (5,6,27,31,32,52,57,94,100,160)],
        'dhfr_bulge_mutants':[i-1 for i in (137,153,155)],
        # BACE
        'bace_ticks':[i-1 for i in range(4,389,10)],
        }
    standard_residue_lists['dhfr_watney_network'] = standard_residue_lists['dhfr_WatneyHammes-SchifferJPhysChemBNetwork_shared'] + standard_residue_lists['dhfr_WatneyHammes-SchifferJPhysChemBNetwork_ecolionly'] + [159,]
    standard_residue_lists['dhfr_networks'] = standard_residue_lists['dhfr_AgarwalJPhysChemBNetwork'] +  standard_residue_lists['dhfr_watney_network']
    standard_residue_lists['dhfr_loops'] = standard_residue_lists['dhfr_m20'] + standard_residue_lists['dhfr_gh'] + standard_residue_lists['dhfr_fg'] + standard_residue_lists['dhfr_cd']
    standard_residue_lists['dhfr_helices'] = []
    for i in 'AB AC AE AF'.split(): standard_residue_lists['dhfr_helices'].extend(standard_residue_lists['dhfr_'+i])
    standard_residue_lists['dhfr_rigid_d1'] = []
    for i in 'BB BC BD BE AC AE AF'.split(): standard_residue_lists['dhfr_rigid_d1'].extend(standard_residue_lists['dhfr_'+i])

    def parse_one_list(value):
        result = []
        for k in value.strip().split(','):
            if k in standard_residue_lists:
                result += standard_residue_lists[k]
            else:
                if '-' in k:
                    start,stop = k.split('-')
                    result.extend(range(int(start)-1,int(stop)))
                else:
                    result.append(int(k)-1)
        result = sorted(list(set(result)))
        print "Caring about",len(result),"residues",[i+1 for i in result]
        return result
    try:
        if ':' in value:
            left,right = value.split(':')
            return [parse_one_list(left),parse_one_list(right)]
        else:
            return parse_one_list(value)
    except ValueError:
        raise OptionValueError(
            "option %s is neither a known range of residues or a valid list of integers: %s"%(opt,value))

# Define a new option type, a 1-based list of integers.  Just use
# the guts of the zero-based stuff and add one.

def check_onebasedintlist(option,opt,value):
    result = check_zerobasedintlist(option,opt,value)
    return [i+1 for i in result]

# Define a new option type, a comma-separated list of strings.
def check_strlist(option,opt,value):
    try:
        return [i for i in value.strip().split(',')]
    except ValueError:
        raise OptionValueError(
            "option %s: invalid list of strings: %s"%(opt,value))


class MyOption(Option):
    TYPES=Option.TYPES + ("zerobasedintlist","onebasedintlist","strlist",)
    TYPE_CHECKER=copy(Option.TYPE_CHECKER)
    TYPE_CHECKER["zerobasedintlist"] = check_zerobasedintlist
    TYPE_CHECKER["onebasedintlist"] = check_onebasedintlist
    TYPE_CHECKER["strlist"] = check_strlist

#######################################
#
# Adding options to the Optparse parser
#
#######################################

def add_window_options(parser):
    parser.add_option("--start",dest="start",default=500,type="int",
                      help="Time, in ps, to start the windows.  [default: %default]")
    parser.add_option("--stop",dest="stop",default=10500,type="int",
                      help="Time, in ps, to stop the windows.  [default: %default]")
    parser.add_option("--window-size",dest="windowsize",default=1000,type="int",
                      help="Length, in ps, of window size.  [default: %default]")
    parser.add_option("--window-spacing",dest="windowspacing",default=100,type="int",
                      help="Spacing between windows, in ps.  [default: %default]")

def add_standard_options(parser):
    parser.add_option("--structure-name",dest="structurename",
                      default="1rx1",
                      help="Name of your structure.  E.g. 1RX1 or 1SGZ. [default: %default]")
    parser.add_option("--output-dir",dest="outputdir",
                      default="ptraj_files/",
                      help="Directory where we will put our results.  This should be the same as the directory where we put our ptraj files before, and it should contain the .dat files that ptraj outputs.  [default: %default]  The images will go to <outputdir>/images and the html file will be in <outputdir>")

#######################################
#
# Parsing standard options
#
#######################################

def get_desired(options):
    starts = range(options.start,options.stop-options.windowsize+1,options.windowspacing)
    stops = range(options.start+options.windowsize,options.stop+1,options.windowspacing)
    names = [(starts[i] + options.windowsize/2.0)/1000.0 for i in range(len(starts))]
    # names is the centers of the windows.  start + windowsize/2 converted to NS.

    desired = [('%sps'%i,'%sps'%j,'NS%05.2f'%k) for (i,j,k) in zip(starts,stops,names)]
    return desired

#######################################
#
# Running external programs
#
#######################################

def run(prog,args,verbose=True):
    '''
    wrapper to handle spaces on windows.
    prog is the full path to the program.
    args is a string that we will split up for you.
        or a tuple.  or a list. your call.

    return value is (retval,prog_out)

    e.g.

    (retval,prog_out) = run("/bin/ls","-al /tmp/myusername")
    '''
    try:
        import subprocess,tempfile

        if type(args) == type(''):
            args = tuple(args.split())
        elif type(args) in (type([]),type(())):
            args = tuple(args)
        if os.name in 'nt dos'.split():
            # turns out not to be necessary with subprocess
            #prog = r'"%s"'%prog
            pass
        args = (prog,) + args
        output_file = tempfile.TemporaryFile(mode="w+")
        if verbose:
            print "Running",args
        retcode = subprocess.call(args,stdout=output_file.fileno(),stderr=subprocess.STDOUT)
        output_file.seek(0)
        prog_out = output_file.read()
        output_file.close() #windows doesn't do this automatically
        if verbose:
            print "Results were:"
            print "Return value:",retcode
            print "Output:"
            print prog_out
        return (retcode,prog_out)
    except ImportError:
        # probably python <= 2.4
        if type(args) != type(''):
            args = ' '.join(args)
        cmd = prog + ' ' + args
        if verbose:
            print "Running",cmd
        retcode = os.system(cmd)
        # cannot return prog_out via os.system
        if verbose:
            print "Results were:"
            print "Return value:",retcode
            print "Output:"
            print "\tcould not find subprocess module, so no output reported"
        return (retcode,'')

def read_data(fname):
    """
    Given a filename, try to read in data, paying attention to whatever format it's in.

    We'll try the following in order

     - .numpy
           We assume this is a square numpy array, made with array.tofile()
     - .dat
           A .dat file from ptraj.  This can be read in with scipy.io.read_array()
    """
    from scipy.io import read_array
    import numpy as N
    data = None
    sys.stdout.write("reading "+fname+" ")
    sys.stdout.flush()
    if os.path.isfile(fname+'.numpy'):
        sys.stdout.write("as numpy\n")
        sys.stdout.flush()
        f = file(fname+'.numpy')
        data = N.fromstring(f.read())
        f.close()
        one_side = int(N.sqrt(len(data)))
        if abs(len(data) - one_side*one_side) >= 0.1:
            message = "Only square matrices are supported. %s has len %s."%(fname,one_side)
            sys.stdout.write(message)
            return None
        data.shape = one_side,one_side
    elif os.path.isfile(fname+'.numpy.bz2'):
        sys.stdout.write("as numpy.bz2\n")
        sys.stdout.flush()
        bf = bz2.BZ2File(fname+'.numpy.bz2')
        data = N.fromstring(bf.read())
        bf.close()
        one_side = int(N.sqrt(len(data)))
        if abs(len(data) - one_side*one_side) >= 0.1:
            message = "Only square matrices are supported. %s has len %s."%(fname,one_side)
            sys.stdout.write(message)
            return None
        data.shape = one_side,one_side
            
    elif os.path.isfile(fname):
        sys.stdout.write("as dat\n")
        sys.stdout.flush()
        data = read_array(file(fname))
    elif os.path.isfile(fname+'.bz2'):
        sys.stdout.write("as dat.bz2\n")
        sys.stdout.flush()
        data = read_array(bz2.BZ2File(fname+'.bz2')) 
    else:
        sys.stdout.write("COULD NOT FIND (1)%s\n"%fname)
    return data
    
