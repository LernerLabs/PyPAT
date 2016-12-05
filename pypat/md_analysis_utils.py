#!/usr/bin/env python

"""
Tools for analyzing MD trajectories
"""
from scipy import *

import sys,glob,os,bz2
from .tool_utils import read_data

def setup_ca(reference_fname,indexing=1,skip_resis=[]):
    """
    most of the code in my_md_analysis assumes 1-based indexing
    for this, but make_plot wants 0-based indexing.  in the end,
    i should switch everything over to 0-based, but i'll leave
    this in as a hack for now.
    """
    master_residue_list = []
    reference_coords = []
    ca_atom_id_list = []
    skip_atom_id_list = []
    #
    # We have a few naming conventions to deal with, like
    # 1HD2 vs. HD12.  One thing that seems consistent is that
    # the first non-numeric character in the atom name tells
    # us the element.
    #
    # No need for speed here.
    #
    def is_hydro(atom_name):
        c = atom_name[0]
        if c in '1234567890':
            c = atom_name[1]
        return c in 'H'
    def is_nosp(atom_name):
        c = atom_name[0]
        if c in '1234567890':
            c = atom_name[1]
        return c in 'NOSP'
    def is_mainchain(atom_name,resn):
        return (resn in 'NAP'.split()) or (atom_name in 'CA C N O H HA'.split())
    
        
    hydro_atom_id_list = []
    nonhydro_atom_id_list = []
    nosp_atom_id_list = []
    nonnosp_atom_id_list = []
    mainchain_atom_id_list = []
    sidechain_atom_id_list = []
    mainchain_nonhydro_atom_id_list = []
    sidechain_hydro_atom_id_list = []
    sidechain_nonhydro_atom_id_list = []
    f = open(reference_fname)
    generated_atom_id = 1
    for line in f:
        # This is an annoyingly hardcoded bit. we should use a real PDB parser.
        line = line.replace('NDPH',' NDP ')
        parts = line.split()
        if parts[0] not in ['ATOM','HETATM']:
            continue
        #atom,atom_id,atom_name,resn,chain,resi,x,y,z,occupancy,b,elem_name = parts
        atom,atom_id,atom_name,resn,chain,resi,x,y,z = parts[:9]
        #
        # One little hack to deal with missing chain info.  If it's a number,
        # it's not really the chain.
        #
        try:
            junk = float(chain)
            atom,atom_id,atom_name,resn,resi,x,y,z = parts[:8]
        except ValueError:
            pass
        
        x,y,z = list(map(float,(x,y,z)))
        try:
            resi,atom_id = list(map(int,(resi,atom_id)))
        except ValueError:
            print("Trouble getting resi,atom_id from",resi,atom_id,line.strip())
            raise
        if resi in skip_resis:
            skip_atom_id_list.append(atom_id)
        if is_hydro(atom_name):
            hydro_atom_id_list.append(atom_id)
        else:
            nonhydro_atom_id_list.append(atom_id)

        if is_nosp(atom_name):
            nosp_atom_id_list.append(atom_id)
        else:
            nonnosp_atom_id_list.append(atom_id)

        if is_mainchain(atom_name,resn):
            mainchain_atom_id_list.append(atom_id)
            if not is_hydro(atom_name):
                mainchain_nonhydro_atom_id_list.append(atom_id)
        else:
            sidechain_atom_id_list.append(atom_id)
            if is_hydro(atom_name):
                sidechain_hydro_atom_id_list.append(atom_id)
            else:
                sidechain_nonhydro_atom_id_list.append(atom_id)

        if atom_name != 'CA':
            continue
        #
        # 1RX2 claims to model two conformations for ASP 116.  I only see one
        # in the PDB file, and it's conformation A.
        #
        if resn == 'AASP' and resi == 116 and '1RX2' in reference_fname: resn = 'ASP'
        elif resn == 'AASP': print("Unrecognized residue: AASP")
        master_residue_list.append((resi,resn))
        reference_coords.append((x,y,z))
        ca_atom_id_list.append(atom_id)
    f.close()
    results = master_residue_list,reference_coords,ca_atom_id_list,hydro_atom_id_list,nonhydro_atom_id_list,nosp_atom_id_list,nonnosp_atom_id_list,mainchain_atom_id_list,sidechain_atom_id_list,mainchain_nonhydro_atom_id_list,sidechain_hydro_atom_id_list,sidechain_nonhydro_atom_id_list,skip_atom_id_list

    #
    # Now skip the things in the skip list
    #
    _results = []
    for (_i,thing) in enumerate(results):
        if not thing:
            # empty lists
            _results.append(thing)
            continue
        if isinstance(thing[0], type(())):
            # master_residue_list has (1,'MET'), etc.
            thing = [i for i in thing if i[0] not in skip_atom_id_list]
        elif isinstance(thing[0], type(1)):
            thing = [i for i in thing if i not in skip_atom_id_list]
        else:
            print("wha?",thing)
            a = 1/0
        _results.append(thing)
    results = _results
            

    #
    # And now take care of the indexing (0 or 1) problem
    #
    if indexing == 0:
        print("Fixing indexing")
        for thing in results:
            for i in range(len(thing)):
                if isinstance(thing[0], type(())):
                    # master_residue_list has (1,'MET'), etc.
                    thing[i] = (thing[i][0] - 1,thing[i][1])
                elif isinstance(thing[0], type(1)):
                    thing[i] = thing[i] - 1
                else:
                    print("wha?",thing)
                    a = 1/0
    return results

_parse_count = 0
def parse_ca(fname,master_residue_list):
    global _parse_count
    _parse_count += 1
    if divmod(_parse_count,100)[-1] == 0:
        sys.stdout.write('.')
        sys.stdout.flush()
    residue_list = []
    coords = []
    f = open(fname)
    for line in f:
        parts = line.split()
        if parts[0] not in ['ATOM','HETATM']:
            continue
        atom,atom_id,atom_name,resn,chain,resi,x,y,z,occupancy,b = parts
        if atom_name != 'CA':
            continue
        x,y,z = list(map(float,(x,y,z)))
        resi,atom_id = list(map(int,(resi,atom_id)))
        #
        # In order to compare different structures, we cheat a little bit.
        # Here, we say that all Histidine residues will be called HIS
        #
        if resn in ('HIE','HID','HIP'): resn = 'HIS'
        residue_list.append((resi,resn))
        coords.append((x,y,z))
    f.close()
    if residue_list != master_residue_list:
        #print master_residue_list
        #print residue_list
        print("broken",[(i,residue_list[i],master_residue_list[i]) for i in range(len(residue_list)) if residue_list[i] != master_residue_list[i]])
        sys.exit()
    return coords

def parse_all(fname,master_residue_list):
    global _parse_count
    _parse_count += 1
    if divmod(_parse_count,100)[-1] == 0:
        sys.stdout.write(',')
        sys.stdout.flush()
    coords = []
    residue_list = []
    f = open(fname)
    for line in f:
        parts = line.split()
        if parts[0] not in ['ATOM','HETATM']:
            continue
        atom,atom_id,atom_name,resn,chain,resi,x,y,z,occupancy,b = parts
        x,y,z = list(map(float,(x,y,z)))
        resi,atom_id = list(map(int,(resi,atom_id)))
        #
        # In order to compare different structures, we cheat a little bit.
        # Here, we say that all Histidine residues will be called HIS
        #
        if resn in ('HIE','HID','HIP'): resn = 'HIS'
        if atom_name == 'CA':
            residue_list.append((resi,resn))
        coords.append((x,y,z))
    f.close()
    if residue_list != master_residue_list:
        #print master_residue_list
        #print residue_list
        print("broken",[(i,residue_list[i],master_residue_list[i]) for i in range(len(residue_list)) if residue_list[i] != master_residue_list[i]])
        sys.exit()
    return coords


vds_radii = {'H':1.20,
             'C':1.7,
             'N':1.55,
             'O':1.52,
             'F':1.35,
             'P':1.9,
             'S':1.85,
             'Cl':1.8,
             }
def nap16O7phe31min_dist(s,include_hydro):
    '''
    nap 160 O7N is
ATOM   2499  O7N NAP A 160      31.906  43.045  14.292  0.00  0.00              
    phe31 is these:
ATOM    464  N   PHE A  31      29.869  36.434   7.846  0.00  0.00              
ATOM    465  H   PHE A  31      30.728  36.378   8.374  0.00  0.00              
ATOM    466  CA  PHE A  31      29.490  37.827   7.596  0.00  0.00              
ATOM    467  HA  PHE A  31      28.684  38.052   8.295  0.00  0.00              
ATOM    468  CB  PHE A  31      30.654  38.827   7.849  0.00  0.00              
ATOM    469  HB2 PHE A  31      31.411  38.440   7.165  0.00  0.00              
ATOM    470  HB3 PHE A  31      31.011  38.747   8.876  0.00  0.00              
ATOM    471  CG  PHE A  31      30.331  40.229   7.525  0.00  0.00              
ATOM    472  CD1 PHE A  31      29.332  40.919   8.262  0.00  0.00              
ATOM    473  HD1 PHE A  31      28.856  40.455   9.113  0.00  0.00              
ATOM    474  CE1 PHE A  31      28.974  42.223   7.928  0.00  0.00              
ATOM    475  HE1 PHE A  31      28.310  42.806   8.549  0.00  0.00              
ATOM    476  CZ  PHE A  31      29.667  42.830   6.822  0.00  0.00              
ATOM    477  HZ  PHE A  31      29.388  43.825   6.505  0.00  0.00              
ATOM    478  CE2 PHE A  31      30.522  42.121   5.987  0.00  0.00              
ATOM    479  HE2 PHE A  31      31.146  42.680   5.304  0.00  0.00              
ATOM    480  CD2 PHE A  31      30.872  40.855   6.386  0.00  0.00              
ATOM    481  HD2 PHE A  31      31.673  40.345   5.873  0.00  0.00              
ATOM    482  C   PHE A  31      28.903  37.987   6.204  0.00  0.00              
ATOM    483  O   PHE A  31      27.702  38.399   6.116  0.00  0.00              
    '''
    phe_nonhydro_ids = [i - 1 for i in [464,466,468,471,472,474,476,478,480,482,483]]
    phe_hydro_ids = [i - 1 for i in [465,467,469,470,473,475,477,479,481]]
    nap_id = 2499 - 1
    if include_hydro:
        phes = phe_nonhydro_ids + phe_hydro_ids
    else:
        phes = phe_nonhydro_ids

    dists = []
    for phe in phes:
        dists.append(((s[phe][0] - s[nap_id][0])**2 + (s[phe][1] - s[nap_id][1])**2 + (s[phe][2] - s[nap_id][2])**2)**0.5)
    return min(dists)

def ile50leu28min_dist(s,include_hydro,vdw):
    '''
    ile is residue 50, which corresponds to these:
ATOM    789  N   ILE A  50      36.039  51.305   7.008  0.00  0.00              
ATOM    790  H   ILE A  50      35.197  51.604   7.478  0.00  0.00              
ATOM    791  CA  ILE A  50      35.920  50.306   5.913  0.00  0.00              
ATOM    792  HA  ILE A  50      36.626  49.526   6.197  0.00  0.00              
ATOM    793  CB  ILE A  50      34.548  49.572   5.812  0.00  0.00              
ATOM    794  HB  ILE A  50      34.634  48.734   5.122  0.00  0.00              
ATOM    795  CG2 ILE A  50      34.338  48.759   7.152  0.00  0.00              
ATOM    796 1HG2 ILE A  50      33.464  48.123   7.010  0.00  0.00              
ATOM    797 2HG2 ILE A  50      35.292  48.296   7.399  0.00  0.00              
ATOM    798 3HG2 ILE A  50      33.946  49.490   7.859  0.00  0.00              
ATOM    799  CG1 ILE A  50      33.333  50.451   5.455  0.00  0.00              
ATOM    800 2HG1 ILE A  50      33.017  51.019   6.330  0.00  0.00              
ATOM    801 3HG1 ILE A  50      33.693  51.294   4.863  0.00  0.00              
ATOM    802  CD1 ILE A  50      32.267  49.782   4.626  0.00  0.00              
ATOM    803 1HD1 ILE A  50      31.563  50.542   4.285  0.00  0.00              
ATOM    804 2HD1 ILE A  50      31.713  49.138   5.308  0.00  0.00              
ATOM    805 3HD1 ILE A  50      32.795  49.237   3.843  0.00  0.00              
ATOM    806  C   ILE A  50      36.481  50.726   4.555  0.00  0.00              
ATOM    807  O   ILE A  50      36.965  49.871   3.850  0.00  0.00
    leu28 is these:
ATOM    411  N   LEU A  28      34.634  36.196   7.999  0.00  0.00              
ATOM    412  H   LEU A  28      35.557  35.797   8.082  0.00  0.00              
ATOM    413  CA  LEU A  28      34.258  36.995   6.836  0.00  0.00              
ATOM    414  HA  LEU A  28      33.750  37.861   7.263  0.00  0.00              
ATOM    415  CB  LEU A  28      35.554  37.350   6.057  0.00  0.00              
ATOM    416  HB2 LEU A  28      35.988  36.377   5.826  0.00  0.00              
ATOM    417  HB3 LEU A  28      35.319  37.721   5.059  0.00  0.00              
ATOM    418  CG  LEU A  28      36.585  38.145   6.693  0.00  0.00              
ATOM    419  HG  LEU A  28      36.899  37.709   7.641  0.00  0.00              
ATOM    420  CD1 LEU A  28      37.845  38.434   5.827  0.00  0.00              
ATOM    421 1HD1 LEU A  28      38.434  39.296   6.138  0.00  0.00              
ATOM    422 2HD1 LEU A  28      38.500  37.569   5.925  0.00  0.00              
ATOM    423 3HD1 LEU A  28      37.485  38.637   4.819  0.00  0.00              
ATOM    424  CD2 LEU A  28      36.097  39.603   6.905  0.00  0.00              
ATOM    425 1HD2 LEU A  28      35.930  40.049   5.925  0.00  0.00              
ATOM    426 2HD2 LEU A  28      36.688  40.264   7.539  0.00  0.00              
ATOM    427 3HD2 LEU A  28      35.166  39.533   7.469  0.00  0.00              
ATOM    428  C   LEU A  28      33.324  36.222   5.961  0.00  0.00              
ATOM    429  O   LEU A  28      32.443  36.867   5.383  0.00  0.00              
    '''
    ile_nonhydro_ids = [i - 1 for i in [789,791,793,795,799,802,806,807]]
    ile_hydro_ids = [i - 1 for i in [790,792,794,796,797,798,800,801,803,804,805]]
    leu_nonhydro_ids = [i - 1 for i in [411,413,415,418,420,424,428,429]]
    leu_hydro_ids = [i - 1 for i in [412,414,416,417,419,421,422,423,425,426,427]]
    if include_hydro:
        iles = ile_nonhydro_ids + ile_hydro_ids
        leus = leu_nonhydro_ids + leu_hydro_ids
    else:
        iles = ile_nonhydro_ids
        leus = leu_nonhydro_ids

    dists = []
    for ile in iles:
        for leu in leus:
            dists.append(((s[ile][0] - s[leu][0])**2 + (s[ile][1] - s[leu][1])**2 + (s[ile][2] - s[leu][2])**2)**0.5)
    return min(dists)
    
    
def ile50leu28ca_dist(s):
    return ((s[49][0] - s[27][0])**2 + (s[49][1] - s[27][1])**2 + (s[49][2] - s[27][2])**2)**0.5

def ile50asp27ca_dist(s):
    return ((s[49][0] - s[26][0])**2 + (s[49][1] - s[26][1])**2 + (s[49][2] - s[26][2])**2)**0.5

def bynumber_dist(s,n1,n2):
    return ((s[n1][0] - s[n2][0])**2 + (s[n1][1] - s[n2][1])**2 + (s[n1][2] - s[n2][2])**2)**0.5

def simple_rmsd(s1,s2):
    """
    rmsd = sqrt((1/N)*sum((x_n - y_n)^2))
    """
    N = len(s1)
    _range = list(range(N))
    deviations = [(s1[i][0] - s2[i][0])**2 +
                  (s1[i][1] - s2[i][1])**2 +
                  (s1[i][2] - s2[i][2])**2 for i in _range]
    return ((1/float(N))*sum(deviations))**0.5

def write_rmsds_and_ile50leu28ca_dists(ref_fname,
                                       target_pattern,
                                       rmsd_out_fname,
                                       ileleu_out_fname,
                                       ileleu_min_out_fname,
                                       ileleu_min_hydro_out_fname,
                                       ileleu_min_vdw_out_fname,
                                       ileleu_min_vdw_hydro_out_fname,
                                       nap160phe31_out_fname,
                                       nap160phe31_min_out_fname,
                                       nap160phe31_min_hydro_out_fname,
                                       ileasp_out_fname,
                                       out_dir):
    '''
    Wrapper function to write out all of the RMSDS and distances I care about

    target_pattern: we glob this together to get the target files
    '''
    structure_fnames = glob.glob(target_pattern)
    if not structure_fnames:
        print("I found zero files matching this pattern:",target_pattern)
    master_residue_list,reference_coords,ca_atom_id_list,hydro_atom_id_list,nonhydro_atom_id_list,nosp_atom_id_list,nonnosp_atom_id_list,mainchain_atom_id_list,sidechain_atom_id_list,mainchain_nonhydro_atom_id_list,sidechain_hydro_atom_id_list,sidechain_nonhydro_atom_id_list,skip_atom_id_list = setup_ca(ref_fname)
    every_atom_coords = [parse_all(fname,master_residue_list) for fname in structure_fnames]
    all_coords = [parse_ca(fname,master_residue_list) for fname in structure_fnames]
    if 1:
        # individual distances
        #
        # ILE50 LEU20 CA dists don't care about the loops, so do them first
        #


        ile50leu28min_dists = [ile50leu28min_dist(coords,include_hydro=False,vdw=False) for coords in every_atom_coords]
        f = open(os.path.join(out_dir,ileleu_min_out_fname),'w')
        for ileleud in ile50leu28min_dists: f.write('%s\n'%ileleud)
        f.close()

        ile50leu28min_dists = [ile50leu28min_dist(coords,include_hydro=True,vdw=False) for coords in every_atom_coords]
        f = open(os.path.join(out_dir,ileleu_min_hydro_out_fname),'w')
        for ileleud in ile50leu28min_dists: f.write('%s\n'%ileleud)
        f.close()

        ile50leu28min_dists = [ile50leu28min_dist(coords,include_hydro=False,vdw=True) for coords in every_atom_coords]
        f = open(os.path.join(out_dir,ileleu_min_vdw_out_fname),'w')
        for ileleud in ile50leu28min_dists: f.write('%s\n'%ileleud)
        f.close()

        ile50leu28min_dists = [ile50leu28min_dist(coords,include_hydro=True,vdw=True) for coords in every_atom_coords]
        f = open(os.path.join(out_dir,ileleu_min_vdw_hydro_out_fname),'w')
        for ileleud in ile50leu28min_dists: f.write('%s\n'%ileleud)
        f.close()

        nap160O7phe31cg_dists = [bynumber_dist(coords,2499-1,471-1) for coords in every_atom_coords]
        f = open(os.path.join(out_dir,nap160phe31_out_fname),'w')
        for npd in nap160O7phe31cg_dists: f.write('%s\n'%npd)
        f.close()

        nap160O7phe31_min_dists = [nap16O7phe31min_dist(coords,include_hydro=False) for coords in every_atom_coords]
        f = open(os.path.join(out_dir,nap160phe31_min_out_fname),'w')
        for npd in nap160O7phe31_min_dists: f.write('%s\n'%npd)
        f.close()

        nap160O7phe31_min_hydro_dists = [nap16O7phe31min_dist(coords,include_hydro=True) for coords in every_atom_coords]
        f = open(os.path.join(out_dir,nap160phe31_min_hydro_out_fname),'w')
        for npd in nap160O7phe31_min_hydro_dists: f.write('%s\n'%npd)
        f.close()
        
        if 0:

            ile50leu28ca_dists = [ile50leu28ca_dist(coords) for coords in all_coords]
            f = open(os.path.join(out_dir,ileleu_out_fname),'w')
            for ileleud in ile50leu28ca_dists: f.write('%s\n'%ileleud)
            f.close()


            ile50asp27ca_dists = [ile50asp27ca_dist(coords) for coords in all_coords]
            f = open(os.path.join(out_dir,ileasp_out_fname),'w')
            for ileaspd in ile50asp27ca_dists: f.write('%s\n'%ileaspd)
            f.close()

    if 0:

        #
        # Now do the loops
        #
        standard_parts = 'm20 all FG CD GH noloops subdomain1 subdomain2 substrate_binding substrate_binding2'.split() 
        residue_chunks = '0-20 20-40 40-60 60-80 80-100 100-120 120-140 140-160'.split()
        anti_correlated_parts = 'm20 40-50 59-68 71-74 95-97'.split() + '116-124 142-149'.split() # 116-124 is correlated with 40-50 and 59-68, 142-149 is correlated with 40-50
        correlated_parts = '110-114 134-141 107-113 149-156 132-142 144-157'.split() + '115-125 1-10 58-62 43-49 38-63 42-51'.split()
        #for structure_parts in correlated_parts + anti_correlated_parts + standard_parts + residue_chunks:
        for structure_parts in 'substrate_binding2'.split():
            this_fname = os.path.join(out_dir,os.path.splitext(rmsd_out_fname)[0] + '_' + structure_parts + os.path.splitext(rmsd_out_fname)[1])
            print("Doing",this_fname)
            structure_parts_map = {# standard parts
                                   'all':(0,162),
                                   'm20':(8,24),
                                   'FG' :(115,132),
                                   'CD' :(63,71),
                                   'GH' :(141,150),
                                   'subdomain2':(37,106),
                                   #Subdomain 2  38-106 (Sawaya/Kraut, Fig 9B)
                                   # residue chunks
                                   '0-20':(0,20),
                                   '20-40':(20,40),
                                   '40-60':(40,60),
                                   '60-80':(60,80),
                                   '80-100':(80,100),
                                   '100-120':(100,120),
                                   '120-140':(120,140),
                                   '140-160':(140,160),
                                   # Anti-correlated parts
                                   '40-50':(39,50),
                                   '59-68':(58,68),
                                   '71-74':(70,74),
                                   '95-97':(94,97),
                                   '116-124':(115,124),
                                   '142-149':(141,149),
                                   # Correlated parts
                                   '110-114':(109,114),
                                   '134-141':(133,141),
                                   '107-113':(106,113),
                                   '149-156':(148,156),
                                   '132-142':(131,142),
                                   '144-157':(143,157),
                                   '115-125':(114,125),
                                   '1-10':(0,10),
                                   '58-62':(57,62),
                                   '43-49':(42,49),
                                   '38-63':(37,63),
                                   '42-51':(41,51),
                                   }
            if structure_parts == 'noloops':
                bad_idxs = list(range(*structure_parts_map['m20'])) + list(range(*structure_parts_map['FG'])) + list(range(*structure_parts_map['CD'])) + list(range(*structure_parts_map['GH']))
                good_idxs = [i for i in range(*structure_parts_map['all']) if i not in bad_idxs]
                print("structure_parts",structure_parts,"residues",good_idxs)
                rmsds = [simple_rmsd([reference_coords[i] for i in good_idxs if i < len(reference_coords)],[coords[i] for i in good_idxs if i < len(coords)]) for coords in all_coords]
            elif structure_parts == 'subdomain1':
                #Subdomain 1  1-37, 107-159 (Sawaya/Kraut, Fig 9B)
                rmsds = [simple_rmsd(reference_coords[0:37]+reference_coords[106:159],coords[0:37]+coords[106:159]) for coords in all_coords]
            elif structure_parts == 'substrate_binding':
                # my best guess at reproducing fig. 1 from Boehr is 27-40,50-59,6-7,112-113
                rmsds = [simple_rmsd(reference_coords[26:40]+reference_coords[49:59]+reference_coords[5:7]+reference_coords[111:113],
                                     coords[26:40]+coords[49:59]+coords[5:7]+coords[111:113]) for coords in all_coords]
            elif structure_parts == 'substrate_binding2':
                # From Osborne2003, 27+36+37+54+96+5+6+26+27+28+29+36+37+50+51+52+54+57
                good_idxs = [i - 1 for i in (5,6,26,27,28,29,36,37,50,51,52,54,57,96)]
                rmsds = [simple_rmsd([reference_coords[i] for i in good_idxs],[coords[i] for i in good_idxs]) for coords in all_coords]
            else:
                print("structure_parts",structure_parts,"residues",structure_parts_map[structure_parts])
                start,stop = structure_parts_map[structure_parts]
                rmsds = [simple_rmsd(reference_coords[start:stop],coords[start:stop]) for coords in all_coords]
            f = open(this_fname,'w')
            for rmsd in rmsds: f.write('%s\n'%rmsd)
            f.close()

def get_resi_to_atom_id_map(pdb_fname,indexing=1):
    """
    most of the code in my_md_analysis assumes 1-based indexing
    for this, but make_plot wants 0-based indexing.  in the end,
    i should switch everything over to 0-based, but i'll leave
    this in as a hack for now.
    """
    atom_id_map = {}
    f = open(pdb_fname)
    generated_atom_id = 1
    for line in f:
        # This is an annoyingly hardcoded bit. we should use a real PDB parser.
        line = line.replace('NDPH',' NDP ')
        parts = line.split()
        if parts[0] not in ['ATOM','HETATM']:
            continue
        #atom,atom_id,atom_name,resn,chain,resi,x,y,z,occupancy,b,elem_name = parts
        atom,atom_id,atom_name,resn,chain,resi,x,y,z = parts[:9]

        #
        # One little hack to deal with missing chain info.  If it's a number,
        # it's not really the chain.
        #
        try:
            junk = float(chain)
            atom,atom_id,atom_name,resn,resi,x,y,z = parts[:8]
        except ValueError:
            pass
        try:
            resi,atom_id = list(map(int,(resi,atom_id)))
        except ValueError:
            print("trouble with resi",resi,"atom_id",atom_id,line.strip())
            raise
        atom_id_map.setdefault(resi,[]).append(atom_id)
    f.close()
    if indexing == 0:
        aim = {}
        for k,v in atom_id_map.items():
            aim[k-1] = [i - 1 for i in v]
        atom_id_map = aim
    return atom_id_map

def get_max_min_ca_resi_data(all_atom_data,atom_id_map,ca_atom_id_list,non_ca_resis=[160,]):
    """
    non_ca_resis is a list of resis that do not contain alpha carbons.
    
    returns everything as a dict.  you may run out of memory here.
    """
    abs_result = zeros((len(atom_id_map),len(atom_id_map)))
    max_result = zeros((len(atom_id_map),len(atom_id_map)))
    min_result = zeros((len(atom_id_map),len(atom_id_map)))
    ca_result  = zeros((len(atom_id_map),len(atom_id_map)))
    #nonhydro_result = zeros((len(nonhydro_atom_id_list),len(nonhydro_atom_id_list)))
    #
    # remember that atom_id and resi are 1-based, whie the array is 0-based.
    # practically, that means that we use 1-based indices when looking things
    # up in atom_id_map and 0-based indices when looking things up in
    # all_atom_data and XXX_results.
    #
    for r1 in atom_id_map:
        for r2 in atom_id_map:
            atom_ids1 = atom_id_map[r1]
            atom_ids2 = atom_id_map[r2]
            try:
                all_vals = [all_atom_data[i1-1][i2-1] for i1 in atom_ids1 for i2 in atom_ids2]
            except IndexError:
                print("i1",i1,"i2",i2)
                print("all_atom_data.keys",all_atom_data.shape)
                print("all_atom_data[i1-1].keys",all_atom_data[i1-1].shape)
                print("r1,r2",r1,r2)
                print("atom_ids1",atom_ids1)
                print("atom_ids2",atom_ids2)
                raise

            if (r1 in non_ca_resis) or (r2 in non_ca_resis):
                #
                # Special case for NAP.  That's resi 160 and has no CA.
                # In general, we can probably do r1-1 < len(ca_atom_list), etc.
                #
                continue
            else:
                ca1 = ca_atom_id_list[r1-1]
                ca2 = ca_atom_id_list[r2-1]
                ca_result[r1-1][r2-1] = all_atom_data[ca1-1][ca2-1]

            ma = max(all_vals)
            mi = min(all_vals)
            max_result[r1-1][r2-1] = ma
            min_result[r1-1][r2-1] = mi
            if abs(ma) > abs(mi):
                abs_result[r1-1][r2-1] = ma
            else:
                abs_result[r1-1][r2-1] = mi
    return {'abs':abs_result,
            'max':max_result,
            'min':min_result,
            'ca':  ca_result,
            }
def write_max_min_ca_resi_versions(fname,ref_pdb_fname,non_ca_resis=[160,],overwrite=False):
    '''
    fname is the .dat file, a simple matrix file output by ptraj.
    ref_pdb_fname is used to get atom names and resi names, so that we can
                  properly select the max,min,abs,ca atoms from each resi.
    overwrite tells us whether or not to overwrite existing files.
    non_ca_resis is a list of resis that do not contain alpha carbons.
    '''
    fnames = {'max':'_resi_max'.join(os.path.splitext(fname))+'.bz2',
              'min':'_resi_min'.join(os.path.splitext(fname))+'.bz2',
              'abs':'_resi_abs'.join(os.path.splitext(fname))+'.bz2',
              'ca' :'_resi_ca'.join( os.path.splitext(fname))+'.bz2',
              }
    pre_existing_files = False
    for fn in list(fnames.values()):
        if os.path.isfile(fn):
            if overwrite:
                print(fn,"already exists, will overwrite")
            else:
                print(fn,"already exists, will not overwrite")
                pre_existing_files = True
    if pre_existing_files and not overwrite:
        return
    atom_id_map = get_resi_to_atom_id_map(ref_pdb_fname)
    master_residue_list,reference_coords,ca_atom_id_list,hydro_atom_id_list,nonhydro_atom_id_list,nosp_atom_id_list,nonnosp_atom_id_list,mainchain_atom_id_list,sidechain_atom_id_list,mainchain_nonhydro_atom_id_list,sidechain_hydro_atom_id_list,sidechain_nonhydro_atom_id_list,skip_atom_id_list = setup_ca(ref_pdb_fname)
    sys.stdout.write('  reading in data\n')
    sys.stdout.flush()
    all_atom_data = read_data(fname)
         
    sys.stdout.write('  ready to calculate\n')
    sys.stdout.flush()
    data = get_max_min_ca_resi_data(all_atom_data,atom_id_map,ca_atom_id_list,non_ca_resis=non_ca_resis)
    for d in data:
        if fnames[d].endswith('.bz2'):
            f = bz2.BZ2File(fnames[d],'w')
        else:
            f = open(fnames[d],'w')
            
        print("writing",fnames[d])
        io.write_array(f,data[d])
        f.close()
    
def ThrLett_to_OneLett(resi, suppress_alert = True):
    """
    Usage:  ThrLett_to_OneLett(resi, suppress_alert = True)

    This function receives resi, the three-letter code 
    for an amino acid and returns its one-letter code.
    If the three-letter code is not recognized, it is 
    simply returned.  If suppress_alert is False, 
    a message will be printed alerting the user that 
    the residue was not recognized
    """

    resiu = resi.upper()
    if resiu == 'ALA':
        return 'A'
    elif resiu == 'ARG' or resiu == 'ARN':
        return 'R'
    elif resiu == 'ASN':
        return 'N'
    elif resiu == 'ASP' or resiu == 'ASH':
        return 'D'
    elif resiu == 'CYS' or resiu == 'CYX' or \
         resiu == 'CYM':
        return 'C'
    elif resiu == 'GLN':
        return 'Q'
    elif resiu == 'GLU' or resiu == 'GLH':
        return 'E'
    elif resiu == 'GLY':
        return 'G'
    elif resiu == 'HIS' or resiu == 'HIE' or \
         resiu == 'HID' or resiu == 'HIP':
        return 'H'
    elif resiu == 'ILE':
        return 'I'
    elif resiu == 'LEU':
        return 'L'
    elif resiu == 'LYS' or resiu == 'LYN':
        return 'K'
    elif resiu == 'MET':
        return 'M'
    elif resiu == 'PHE':
        return 'F'
    elif resiu == 'PRO':
        return 'P'
    elif resiu == 'SER':
        return 'S'
    elif resiu == 'THR':
        return 'T'
    elif resiu == 'TRP':
        return 'W'
    elif resiu == 'TYR':
        return 'Y'
    elif resiu == 'VAL':
        return 'V'
    else:
        if not suppress_alert:
            print(resi, "not recognized as residue.  Returning", resi)
        return resi 


def OneLett_to_ThrLett(resi, cap = 'standard', suppress_alert = True):
    """
    Usage:  OneLett_to_ThrLett(resi, cap = 'standard', 
               suppress_alert = True)

    This function receives resi, a one-letter amino acid
    code, and returns the three letter amino acid code.
    If cap(italization) is standard, then only the first
    letter will be capitalized.  If cap is all, then all 
    letters will be capitalized (as in a PDB file).  If
    the one-letter amino acid is not recognized, it is 
    simply returned.  If option suppress_error is False,
    a message will also be printed to alert to user that
    the code was not recognized.
    """
    if resi == 'A':
        res3 = 'Ala'
    elif resi == 'C':
        res3 = 'Cyx'
    elif resi == 'D':
        res3 = 'Asp'
    elif resi == 'E':
        res3 = 'Glu'
    elif resi == 'F':
        res3 = 'Phe'
    elif resi == 'G':
        res3 = 'Gly'
    elif resi == 'H':
        res3 = 'His'
    elif resi == 'I':
        res3 = 'Ile'
    elif resi == 'K':
        res3 = 'Lys'
    elif resi == 'L':
        res3 = 'Leu'
    elif resi == 'M':
        res3 = 'Met'
    elif resi == 'N':
        res3 = 'Asn'
    elif resi == 'P':
        res3 = 'Pro'
    elif resi == 'Q':
        res3 = 'Gln'
    elif resi == 'R':
        res3 = 'Arg'
    elif resi == 'S':
        res3 = 'Ser'
    elif resi == 'T':
        res3 = 'Thr'
    elif resi == 'V':
        res3 = 'Val'
    elif resi == 'W':
        res3 = 'Trp'
    elif resi == 'Y':
        res3 = 'Tyr'
    else:
        if not suppress_alert:
            print(resi, "not recognized as residue.  Returning", resi)
        res3 = resi

    if cap == 'all':
        res3 = res3.upper()
    elif cap != 'standard':
        print('Ignoring invalid option for cap:', cap)

    return res3


def get_resinum_to_resi_map(resiname_file, offset = 0, indexing = 1, aa_code = 3):
    """
    This function returns a dictionary that relates the residue
    number for a given system to a residue name.  A PDB or prmtop file 
    for the system of interest is required.  The string that comprises
    the residue name (the values of the dictionary) is the amino
    acid code (in three-letter or one-letter format; default three-
    letter) followed by the residue number (e.g., Thr72 or T72).

    An optional offset can be passed to the function.  This 
    changes the residue number in the string, which is useful if 
    the numbering in the PDB file is not the preferred numbering 
    for the system.  

    For most indexed applications, the first residue is
    numbered as 0 instead of 1.  A zero-based dictionary can be
    set by passing indexing = 0 to the function.  The default
    is 1.  

    If the PDB file input to the function is not
    found, a message will alert the user.  In this case, the
    values in the dictionary are simply the numbers from 1 to 
    9999 (the maximum residue number of a PDB file).  
    """
    resi_map = {}

    if resiname_file == None:
        print('Warning:  No prmtop or PDB file given.\n' + \
              '  No residue number information will be presented.')
    for i in range(10000):
        resi_map[i] = str(i)
    return resi_map

    try:
        f = file(resiname_file)
    except IOError:
        print('Warning:  Could not open ' + resiname_file + '.\n' + \
              '  No residue number information will be presented.')
    for i in range(10000):
        resi_map[i] = str(i)
    return resi_map

    # If the file is a prmtop file...

    if not resiname_file.endswith('.pdb'):
        resi_num = 1
        
        residue_section = False
        for line in f:
            if line.startswith('%FLAG RESIDUE_POINTER'):
                break
            if line.startswith('%FLAG RESIDUE_LABEL'):
                residue_section = True
            if not residue_section or line.startswith('%F'):
                continue
            else:
                residue_names = line.split()
        for resi_name in residue_names:
            if aa_code == 1:
                resi_name = ThrLett_to_OneLett(resi_name)
                resi_name = resi_name.capitalize() + str(resi_num + offset)
                resi_map[resi_num + indexing - 1] = resi_name
                resi_num += 1

    # If the file is a PDB file...

    else:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            resi_name = line[17:21].strip()
            resi_num = int(line[22:26].strip())
            if aa_code == 1:
                resi_name = ThrLett_to_OneLett(resi_name)
                resi_name = resi_name.capitalize() + str(resi_num + offset)
                resi_map[resi_num + indexing - 1] = resi_name
    
    f.close()

    if not resi_map:
        print("Warning: Could not extract residue information from prmtop or PDB file.\n")
        print("  No residue number information will be presented.")
    for i in range(10000):
        resi_map[i] = str(i)
    return resi_map
    
    return resi_map

