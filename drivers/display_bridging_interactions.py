#!/cluster/home2/mglerner/anaconda3/bin/python
#!/usr/bin/env python

import sys,os

if __name__ == '__main__':
    from pypat.hbond import pymol_hbond_analysis
    from pypat import tool_utils

    from optparse import OptionParser    

    parser = OptionParser(option_class=tool_utils.MyOption)

    parser.add_option("--name",'-n',dest="name",default="nrna",
                      help="Trajectory and topology must be named name.trj and name.top respectively. [default: %default]"
                      )
    parser.add_option("--timestep",'-t',dest='timestep',default=5,
                      help="trajectory timestep in picoseconds. [default: %default]",
                      type='int',
                      )
    parser.add_option('--min-dwell-time','-m',dest='minrequireddwelltime',default=3,
                      help='minimum required dwell time. [default: %default]',
                      type='int',
                      )
    parser.add_option('--looseness','-l',dest='looseness',default=2,
                      help='looseness. [default: %default]',
                      type='int',
                      )
    parser.add_option("--dist-cutoff",'-d',dest="distcutoff",default=3.5,
                      help="Heavy atom to heavy atom distance cutoff. [default: %default]",
                      type='float',
                      )
    parser.add_option('--angle-cutoff','-a',dest='anglecutoff',default=0.0,
                      help="Angle cutoff. If heavy:hydro:heavy angle must be greater than this. [default: %default]",
                      type='float',
                      )
    parser.add_option('--min-occ','-o',dest='minocc',default=0.4,
                      help="Minimum percentage of the trajectory for which this interaction must be occupied. [default: %default]",
                      type='float',
                      )
    parser.add_option('--dir','-r',dest='dir',default='.',
                      help="Directory in which the name_hbond_*_*.txt files are located. [default: %default]",
                      )
    parser.add_option('--resi-criteria','-R',dest='resi_criteria',default=None,
                      help="Restrict the output to BWIs where at least one side involves a residue in this list. [default: %default]",
                      type='onebasedintlist',
                      )

    options,args = parser.parse_args()

    


    import glob
    d = options.dir
    struct = options.name
    fnames = glob.glob(os.path.join(d,'%s_hbond_?_?.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_?_??.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_??_??.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_??_???.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_???_???.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_???_????.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_????_????.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_?????_?????.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_??????_??????.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_???????_???????.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_????????_????????.txt'%struct))
    fnames += glob.glob(os.path.join(d,'%s_hbond_?????????_?????????.txt'%struct))



    a = pymol_hbond_analysis.get_hbond_trajectories(structure=options.name,
                                                    timestep=options.timestep,
                                                    fnames=fnames,
                                                    combination_method='loose',
                                                    min_required_dwell_time=options.minrequireddwelltime,
                                                    looseness=options.looseness,
                                                    dist_cutoff=options.distcutoff,
                                                    angle_cutoff=options.anglecutoff,
                                                    ) # skip the last one because we're probably writing to it right now
    #a = get_hbond_trajectories('1rx1',5,fnames1[:50]) # skip the last one because we're probably writing to it right now

    #nap_including_resis = list((160,)) # NAP
    #m20_including_resis = list(range(9,25)) # m20
    #tunnel_including_resis = list((5,7,8,30,33,34,92,111,112,113,114,137,153,155)) # tunnel
    #tunnel_surface_including_resis =  list((30,33,111,137,153,155)) # tunnel surface
    #no_including_resis = None
    #including_resis = tunnel_including_resis + nap_including_resis
    #including_resis = no_including_resis
    print(a.get_trajectory_string(#minocc=0.05,
                                  #minocc=0.60,
                                  minocc=options.minocc,
                                  numchunks=50,
                                  including_resis=options.resi_criteria,
                                  #sort_by='average_dwell_time',
                                  sort_by='occ',
                                  ))

    
