#!/usr/bin/env python

import sys


if __name__ == '__main__':
    from pypat.hbond import pymol_hbond_analysis
    from optparse import OptionParser    

    usage = """
    Run this like:

    pymol -qcr collect_water_bridges.py -- --name=myprotein --dist-cutoff=3.5

    Do not forget the double dashes after the script name.
    """
    parser = OptionParser(usage=usage)

    parser.add_option("--name",'-n',dest="name",default="nrna",
                      help="Trajectory and topology must be named name.trj and name.top respectively. [default: %default]"
                      )
    parser.add_option("--chunk-size",'-c',dest="chunksize",default=10,
                      help="How many MD steps to process at a time. If you do too many at a time, PyMOL will slow down. Too few, and you're wasting time starting/stopping PyMOL. [default: %default]",
                      type='int',
                      )
    parser.add_option("--num-steps",'-s',dest="numsteps",default=50,
                      help="number of steps in your MD trajectory. [default: %default]",
                      type='int',
                      ##FIXME: can be extracted from PyMOL
                      )
    parser.add_option("--dist-cutoff",'-d',dest="distcutoff",default=4.0,
                      help="Heavy atom to heavy atom distance cutoff. [default: %default]",
                      type='float',
                      )
    parser.add_option('--angle-cutoff','-a',dest='anglecutoff',default=0.0,
                      help="Angle cutoff. If heavy:hydro:heavy angle must be greater than this. [default: %default]",
                      type='float',
                      )

    argv = sys.argv
    if '--' in argv:
        argv = argv[argv.index('--')+1:]
    print(argv)
    options,args = parser.parse_args(argv)

    r = list(range(1,options.numsteps+1,options.chunksize))
    if r[-1] != options.numsteps:
        r.append(options.numsteps)
    starts_and_stops = list(zip(r[:-1],r[1:]))
    print(starts_and_stops)
    for (start,stop) in starts_and_stops:
        print()
        print('DOING',start,stop)
        print()
        pymol_hbond_analysis.find_bridging_waters_in_trajectory(name=options.name,
                                                                start=start,
                                                                stop=stop,
                                                                hbond_dist_cutoff=options.distcutoff,
                                                                hbond_angle_cutoff=options.anglecutoff,
                                                                )
