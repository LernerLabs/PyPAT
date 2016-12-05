#!/usr/bin/env python
'''
runs ptraj lots of times.
'''

import os,sys

def run_ptraj_hbond_analysis(desired_params):
    for (start,stop,name) in desired_params:
        for struct in '1ra1'.split():
            print(struct,name)
            ptraj = 'ptraj'
            prmtop = '/users/mlerner/work/src/Dynamics-DHFR/MD_Files/%s/%s_noha22_solv.prmtop.resized'%(struct.upper(),struct)
            ptraj_input_file_dir = '/users/mlerner/work/src/Dynamics-DHFR/MD_Files/%s/TrajectoriesForIndividualNanoseconds/PtrajFiles/'%(struct.upper())
            fname = 'ptraj_hbond_%s_%s.ptraj'%(struct,name)
            ptraj_input_filename =   os.path.join(ptraj_input_file_dir,fname)
            hbond_output_file_dir = '/users/mlerner/work/src/Dynamics-DHFR/MD_Files/%s/TrajectoriesForIndividualNanoseconds/PtrajFiles/'%(struct.upper())
            fname = 'hbond_%s_%s.out'%(struct,name)
            hbond_output_filename = os.path.join(hbond_output_file_dir,fname)
            # Run ptraj
            print("Running ptraj")
            cmd = 'ptraj prmtop < input > output'
            cmd = '%s %s < %s > %s'%(ptraj,prmtop,ptraj_input_filename,hbond_output_filename)
            print(cmd)
            os.system(cmd)

if __name__ == '__main__':

    starts = list(range(500,9501,100))
    stops  = list(range(1500,10501,100))
    names  = [1 + i/10.0 for i in range(len(starts))]
    names = [i for i in names if i > 5.9]
    desired_params = [('%sps'%i,'%sps'%j,'NS%3.1f'%k) for (i,j,k) in zip(starts,stops,names)]
    # 1RA1 needs to start on 4.5
    # 1RX1 needs to start on 8.2

    run_ptraj_hbond_analysis(desired_params)

