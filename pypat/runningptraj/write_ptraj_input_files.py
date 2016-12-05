#!/usr/bin/env python

import os

def write_ptraj_input_files(dir_containing_ptraj_files,
                            ptraj_output_dir,filename_prefix,desired,ps_per_frame,ptraj_header,fname_template,
                            write_covar,
                            write_out_pdb_file=None,pdb_ptraj_header=None):
    """
    dir_containing_ptraj_files is the place where we'll write out the ptraj
                               input files

    ptraj_output_dir is the directory where ptraj will spit the results out

    filename_prefix is the start of the output filenames.  E.g. if you set it
                    to 1ra1, you will get files like
                    <output_dir>/1ra1_1NS_all_atom_covarmat.dat

    desired is a list of starts, stops and names, zipped together.  See the
            examples below.

    ps_per_frame is the number of ps per frame in the trajectory file.

    ptraj_header is the header of the ptraj file.  it should load in the
                 trajectories, strip out all of the atoms you don't care
                 about and align things.  See an example below.

    write_covar if this is true, we'll write out the covariance matrix.  it's
                a big matrix that we don't always use, so we usually don't
                write it out.

    fname_template is ...

    if write_out_pdb_file is not None, we'll write out a PDB file that
                          corresponds to this ptraj run.

    pdb_ptraj_header is the ptraj header we will use to write out that pdb file.
                     It really only needs to read in the first state of the
                     first mdcrd file.
     
    """
    #
    # each frame is ps_per_frame ps, so we can turn the human-readable things in the 'desired'
    # list into the appropriate frames for the ptraj files.
    #
    for (start,stop,name) in desired:
        start = int(start[:-2])
        stop =  int(stop[:-2])
        if start != 1:
            start = int(start/ps_per_frame)
        stop = int(stop/ps_per_frame)
        fname = fname_template % name
        #print fname
        if write_covar:
            txt = ptraj_header +  '''
matrix correl out %s/%s_%s_all_atom_correlmat.dat start %s stop %s byatom
matrix covar out %s/%s_%s_all_atom_covarmat.dat start %s stop %s
matrix correl out %s/%s_%s_byres_correlmat.dat start %s stop %s byres
average  %s/%s_%s_average.pdb start %s stop %s pdb
go
    '''%(ptraj_output_dir,filename_prefix,name,start,stop,
         ptraj_output_dir,filename_prefix,name,start,stop,
         ptraj_output_dir,filename_prefix,name,start,stop,
         ptraj_output_dir,filename_prefix,name,start,stop,
         )
        else:
            txt = ptraj_header +  '''
matrix correl out %s/%s_%s_all_atom_correlmat.dat start %s stop %s byatom
matrix correl out %s/%s_%s_byres_correlmat.dat start %s stop %s byres
average  %s/%s_%s_average.pdb start %s stop %s pdb
go
    '''%(ptraj_output_dir,filename_prefix,name,start,stop,
         ptraj_output_dir,filename_prefix,name,start,stop,
         ptraj_output_dir,filename_prefix,name,start,stop,
         )
            
        file(os.path.join(dir_containing_ptraj_files,fname),'w').write(txt)

    if write_out_pdb_file:
        txt = pdb_ptraj_header + '''
trajout %s pdb
go
'''%write_out_pdb_file
        file(os.path.join(dir_containing_ptraj_files,'write_%s_ref_pdb.ptraj'%filename_prefix),'w').write(txt)
        
        

