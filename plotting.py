#!/usr/bin/env python
"""

TODO:

1) Get rid of unused functions
2) Possibly move the cmap lookup to tool_utils.

"""
from __future__ import division
from sentinel_map import SentinelMap,SentinelNorm
import glob
import matplotlib.numerix.ma as ma
import md_analysis_utils
import sys
import pylab,matplotlib,os
from matplotlib import mlab
import numpy as N
import scipy.io
import pprint
import bz2
import copy
from tool_utils import read_data
cdict = {'red':  (((-1.0+1)/2,32/255, 32/255),
                  ((-0.6+1)/2,32/255, 32/255),
                  ((-0.3+1)/2,32/255, 32/255),
                  ((-0.2+1)/2, 0/255,  0/255),
                  ((-0.1+1)/2, 0/255,  0/255),
                  (( 0.0+1)/2, 0/255,  0/255),
                  (( 0.1+1)/2, 0/255,  0/255),
                  (( 0.2+1)/2, 76/255, 76/255),
                  (( 0.3+1)/2,255/255,255/255),
                  (( 0.4+1)/2,255/255,255/255),
                  (( 0.5+1)/2,255/255,255/255),
                  (( 0.6+1)/2,242/255,242/255),
                  (( 0.7+1)/2,236/255,236/255),
                  (( 0.8+1)/2,236/255,236/255),
                  (( 0.9+1)/2,213/255,213/255),
                  (( 1.0+1)/2,108/255,108/255)),
         'green':(((-1.0+1)/2, 31/255, 31/255),
                  ((-0.6+1)/2, 31/255, 31/255),
                  ((-0.3+1)/2, 31/255, 31/255),
                  ((-0.2+1)/2, 86/255, 86/255),
                  ((-0.1+1)/2,131/255,131/255),
                  (( 0.0+1)/2,128/255,128/255),
                  (( 0.1+1)/2,137/255,137/255),
                  (( 0.2+1)/2,184/255,184/255),
                  (( 0.3+1)/2,244/255,244/255),
                  (( 0.4+1)/2,244/255,244/255),
                  (( 0.5+1)/2,172/255,172/255),
                  (( 0.6+1)/2, 52/255, 52/255),
                  (( 0.7+1)/2,  1/255,  1/255),
                  (( 0.8+1)/2,  0/255,  0/255),
                  (( 0.9+1)/2,  0/255,  0/255),
                  (( 1.0+1)/2, 24/255, 24/255)),
         'blue': (((-1.0+1)/2,109/255,109/255),
                  ((-0.6+1)/2,109/255,109/255),
                  ((-0.3+1)/2,109/255,109/255),
                  ((-0.2+1)/2,160/255,160/255),
                  ((-0.1+1)/2,202/255,202/255),
                  (( 0.0+1)/2,200/255,200/255),
                  (( 0.1+1)/2,144/255,144/255),
                  (( 0.2+1)/2, 45/255, 45/255),
                  (( 0.3+1)/2,  0/255,  0/255),
                  (( 0.4+1)/2,  0/255,  0/255),
                  (( 0.5+1)/2,  0/255,  0/255),
                  (( 0.6+1)/2,  0/255,  0/255),
                  (( 0.7+1)/2, 20/255, 20/255),
                  (( 0.8+1)/2, 22/255, 22/255),
                  (( 0.9+1)/2, 27/255, 27/255),
                  (( 1.0+1)/2, 33/255, 33/255)),
         }
scaled_by_half_cdict = {}
for color in cdict:
    bottom = cdict[color][0]
    top = cdict[color][-1]
    scaled_by_half_cdict[color] = [bottom,]
    for value in cdict[color]:
        scaled_by_half_cdict[color].append([0.25+value[0]/2,value[1],value[2]])
    scaled_by_half_cdict[color].append(top)
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
my_scaled_cmap = matplotlib.colors.LinearSegmentedColormap('my_scaled_colormap',scaled_by_half_cdict,256)
##pprint.pprint( cdict['red'])
##pprint.pprint( scaled_by_half_cdict['red'])


sentinel1 = -10
sentinel2 = -20
rgb1 = (0.,0.,0.)
rgb2 = (1.,1.,1.)
#my_sentinelcmap = SentinelMap(my_cmap,sentinels={sentinel:rgb})
#my_sentinelnorm = SentinelNorm(ignore=[sentinel,],vmin=-1.0,vmax=1.0)
sentinel_maps_and_norms = {'Normal':(SentinelMap(my_cmap,sentinels={sentinel1:rgb1,sentinel2:rgb2}),
                                     SentinelNorm(ignore=[sentinel1,sentinel2,],vmin=-1.0,vmax=1.0)
                                     ),
                           'Scaled':(SentinelMap(my_scaled_cmap,sentinels={sentinel1:rgb1,sentinel2:rgb2}),
                                     SentinelNorm(ignore=[sentinel1,sentinel2,],vmin=-1.0,vmax=1.0)
                                     ),
                           }
    
                                    
##                           'jet'   :(SentinelMap(pylab.cm.jet,sentinels={sentinel1:rgb1,sentinel2:rgb2}),
##                                     SentinelNorm(ignore=[sentinel1,sentinel2,],vmin=-1.0,vmax=1.0)
##                                     ),
##                           }
for map in [m for m in pylab.cm.datad.keys() if not m.endswith("_r")]:
    sentinel_maps_and_norms[map] = (SentinelMap(pylab.cm.get_cmap(map),sentinels={sentinel1:rgb1,sentinel2:rgb2}),
                                     SentinelNorm(ignore=[sentinel1,sentinel2,],vmin=-1.0,vmax=1.0)
                                    )

                                     

def plotcontours_in_several_ways(fname='/Users/mglerner/work/Dynamics-DHFR/MD_Files/1RX1/1rx1_byres_correlmat.dat'):
    data = scipy.io.read_array(file(fname))
    contourset = pylab.contourf(data)
    pylab.colorbar()
    pylab.show()
    pylab.pcolor(data)
    pylab.colorbar()
    pylab.show()
    pylab.pcolormesh(data)
    pylab.colorbar()
    pylab.show()
    del data

def plot_colormesh(fname,step=0.01,cmap=my_scaled_cmap):
    data = scipy.io.read_array(file(fname))
    #a = pylab.contourf(data,(-0.7,-0.2,0,0.2,0.6,1.0))
    a = pylab.contourf(data,pylab.arange(-1,1+step,step),cmap=cmap) # fine
    pylab.colorbar()
    pylab.show()
    del data
def plot_diff_colormesh(fn1,fn2):
    d1 = scipy.io.read_array(file(fn1))
    d2 = scipy.io.read_array(file(fn2))
    a = pylab.contourf(d1-d2,pylab.arange(-1,1,0.05),cmap=my_cmap) #coarse
    pylab.colorbar()
    pylab.show()
    

def plot_as_pcolor(fname):
    data = scipy.io.read_array(file(fname))
    a = pylab.pcolor(data,cmap=my_cmap)
    pylab.colorbar()
    pylab.show()
    del data

    


def plotrms(filenames):
    num_plots = len(filenames)
    pylab.figure(1)
    for i,fn in enumerate(filenames):
        pylab.subplot(num_plots,1,i+1)
        data = scipy.io.read_array(fn)
        pylab.plot([d[0] for d in data],[d[1] for d in data])
        pylab.title(fn[:-4])
        del data
    pylab.show()

def plotrms2(struct='1rx1'):
    md_data = scipy.io.read_array('%s_rms_vs_md.txt'%struct)
    ra1_data = scipy.io.read_array('%s_rms_vs_1RA1.txt'%struct)
    rx1_data = scipy.io.read_array('%s_rms_vs_1RX1.txt'%struct)
    rx2_data = scipy.io.read_array('%s_rms_vs_1RX2.txt'%struct)
    rx4_data = scipy.io.read_array('%s_rms_vs_1RX4.txt'%struct)
    rx5_data = scipy.io.read_array('%s_rms_vs_1RX5.txt'%struct)
    rx6_data = scipy.io.read_array('%s_rms_vs_1RX6.txt'%struct)
    #pylab.plot([d[0] for d in md_data],[d[1] for d in md_data],label='%s vs MD'%struct,)
    pylab.plot([d[0] for d in ra1_data],[d[1] for d in ra1_data],label='%s vs 1RA1'%struct,)
    pylab.plot([d[0] for d in rx1_data],[d[1] for d in rx1_data],label='%s vs 1RX1'%struct,)
    pylab.plot([d[0] for d in rx2_data],[d[1] for d in rx2_data],label='%s vs 1RX2'%struct,)
    #pylab.plot([d[0] for d in rx4_data],[d[1] for d in rx4_data],label='%s vs 1RX4'%struct,)
    pylab.plot([d[0] for d in rx5_data],[d[1] for d in rx5_data],label='%s vs 1RX5'%struct,)
    pylab.plot([d[0] for d in rx6_data],[d[1] for d in rx6_data],label='%s vs 1RX6'%struct,)
    pylab.legend()
    pylab.show()

def plot_aligned_rms(struct='1rx1',structure_part='m20',in_dir='/Users/mglerner/work/Dynamics-DHFR/RMSD_Txt_Files/',out_dir = '',include_averages=False,prefix='rmsd_ca',title_mode='long',save_fig=True):
    if 'pymol' in prefix:
        file_title = 'pymol %s vs other structures (aligned, ca, %s)'%(struct,structure_part)
    else:
        file_title = '%s vs other structures (aligned, ca, %s)'%(struct,structure_part)
    if include_averages: file_title = file_title + ' ' + 'ave'
    for s in '1RA1 1RX1 1RX2 1RX4 1RX5 1RX6'.split():
        fname = os.path.join(in_dir,prefix+'_kellyaligned_%s_vs_%s_%s.txt'%(struct,s,structure_part))
        print "Trying to load",fname
        data = scipy.io.read_array(fname)
        if title_mode == 'long':
            title = '%s vs %s (aligned, ca, %s)'%(struct,s,structure_part)
        elif title_mode == 'short':
            title = 'vs. %s'%(s,)
        if include_averages:
            if len(data[0]) == 1:
                ave = sum(data)/len(data)
                title = title + ' Ave. %s'%ave
            elif len(data[0]) == 2:
                ave = sum([i[1] for i in data])/len([i[1] for i in data])
                title = title + ' Ave. %s'%ave
        if len(data.shape) == 1:
            pylab.plot([i*5/1000.0 for i in range(len(data))],data,label=title)
        elif len(data[0]) == 2:
            pylab.plot([i*5/1000.0 for i in range(len(data))],[i[1] for i in data],label=title)
        else:
            sys.exit('what?')
        del data
    pylab.figlegend()
    a = pylab.axis()
    pylab.xticks([i/2.0 for i in range(24)])
    pylab.xticks([0,1,2,3,4,5,6,7,8,9,10],size='larger')
    pylab.yticks(size='larger')
    # the legend goes from 4.5 to 6 on a scale of 0-6, so about 25%
    pylab.axis((0,10.505,a[2],a[3]+0.25*(a[3]-a[2]))) # add one on the y axis so there's room for the legend.
    pylab.grid(b=True)
    pylab.xlabel('Time (ns)',size='larger')
    pylab.rc('text',usetex=True)
    pylab.ylabel(r'RMSD ($\AA{}$)',size='larger')
    if save_fig:
        pylab.savefig(os.path.join('RMSDImages',out_dir,file_title+'.png'))
    pylab.show()
    pylab.clf()
    
def plot_ileleu(struct='1rx1'):
    for s in '1RA1 1RX1 1RX2 1RX4 1RX5 1RX6'.split():
        data = scipy.io.read_array(file('ile50leu20_ca_kellyaligned_%s_vs_%s.txt'%(struct,s)))
        pylab.plot(range(len(data)),data,label='ILE LEU %s vs %s (aligned, ca)'%(struct,s))
        del data
    pylab.legend()
    pylab.show()
    

def plotrms3():
    ra1_data = scipy.io.read_array('rms_vs_1RA1.txt')
    ra1_ca_data = scipy.io.read_array('rms_vs_1RA1_CA.txt')
    rx2_data = scipy.io.read_array('rms_vs_1RX2.txt')
    rx2_ca_data = scipy.io.read_array('rms_vs_1RX2_CA.txt')
    #pylab.plot([d[0] for d in ra1_data],[d[1] for d in ra1_data],label='vs 1RA1')
    pylab.plot([d[0] for d in ra1_ca_data],[d[1] for d in ra1_ca_data],label='vs 1RA1 CA')
    #pylab.plot([d[0] for d in rx2_data],[d[1] for d in rx2_data],label='vs 1RX2')
    pylab.plot([d[0] for d in rx2_ca_data],[d[1] for d in rx2_ca_data],label='vs 1RX2 CA')
    pylab.legend()
    pylab.show()

def OLDplot_ile_leu_dist():
    ra_ave_data = scipy.io.read_array('1ra1ile50leu28ave.dat')
    ra_ca_data = scipy.io.read_array('1ra1ile50leu28ca.dat')
    rx_ave_data = scipy.io.read_array('1rx1ile50leu28ave.dat')
    rx_ca_data = scipy.io.read_array('1rx1ile50leu28ca.dat')
    #pylab.plot([d[0] for d in ra_ave_data],[d[1] for d in ra_ave_data],label='1RA1 Average',)
    #pylab.plot([d[0] for d in ra_ca_data],[d[1] for d in ra_ca_data],label='1RA1 CA',)
    pylab.plot([d[0] for d in rx_ave_data],[d[1] for d in rx_ave_data],label='1RX1 Average',)
    pylab.plot([d[0] for d in rx_ca_data],[d[1] for d in rx_ca_data],label='1RX1 CA',)
    pylab.legend()
    pylab.show()

def make_several_correl_plots(struct,times,cmap,detail_level,save_fig,out_dir='CorrelAndCovarImages/',overwrite=True,
                              dat_dir='/Users/mglerner/work/Dynamics-DHFR/BigDynamicsMoviePtrajFiles/',
                              ref_pdb_fname = '/Users/mglerner/work/Dynamics-DHFR/MD_Files/1RX1/ChainAPDBFilesFromTrajectory/AlignedVs1RX1/1rx1_trajectory_A_snap_0001_trans.pdb',
                              plot_types='ca avg max min abs straight mainheavy allheavy sidechainhbond hbond'.split(),
                              mark_resis=[],
                              highlight=0.0,
                              highlight_mode='positive',
                              skip_resis=[],
                              ticks=True,
                              dpi=200,
                              title=None,
                              ):
    for plot_type in [i for i in 'abs max min avg ca'.split() if i in plot_types]:
        if plot_type == 'avg': fname = os.path.join(dat_dir,'%s/%s_%s_byres_correlmat.dat'%(struct,struct,times))
        else:                  fname = os.path.join(dat_dir,'%s/%s_%s_all_atom_correlmat_resi_%s.dat'%(struct,struct,times,plot_type))
        data = read_data(fname)
        if data is None:
            continue
        _make_one_correl_plot(data,struct,times,cmap,detail_level,plot_type,save_fig,dat_dir,out_dir,overwrite,ref_pdb_fname,mark_resis,highlight,highlight_mode,skip_resis,ticks,dpi,title)
        del data

    longer_plot_types = [i for i in 'straight mainheavy allheavy sidechainhbond hbond'.split() if i in plot_types]
    if longer_plot_types:
        fname = os.path.join(dat_dir,'%s/%s_%s_all_atom_correlmat.dat'%(struct,struct,times))
        data = read_data(fname)
        if data is None:
            return
        for plot_type in longer_plot_types:
            _data = copy.copy(data)
            _make_one_correl_plot(data,struct,times,cmap,detail_level,plot_type,save_fig,dat_dir,out_dir,overwrite,ref_pdb_fname,mark_resis,highlight,highlight_mode,skip_resis,ticks,dpi,title)
            del _data
        del data


def make_one_correl_plot(struct,times,cmap,detail_level,plot_type,save_fig,dat_dir,out_dir,overwrite,ref_pdb_fname,mark_resis=[],highlight=0.0,highlight_mode='positive',skip_resis=[],ticks=True,dpi=200,title=None):
    if plot_type in 'straight mainheavy allheavy sidechainhbond hbond'.split(): fname = os.path.join(dat_dir,'%s/%s_%s_all_atom_correlmat.dat'%(struct,struct,times))
    elif plot_type == 'avg':                                           fname = os.path.join(dat_dir,'%s/%s_%s_byres_correlmat.dat'%(struct,struct,times))
    else:                                                              fname = os.path.join(dat_dir,'%s/%s_%s_all_atom_correlmat_resi_%s.dat'%(struct,struct,times,plot_type))
    data = read_data(fname)
    if data is None: return

    _make_one_correl_plot(data,struct,times,cmap,detail_level,plot_type,save_fig,dat_dir,out_dir,overwrite,ref_pdb_fname,mark_resis,highlight,highlight_mode,skip_resis,ticks,dpi,title)

def _make_one_correl_plot(data,struct,times,cmap,detail_level,plot_type,save_fig,dat_dir,out_dir,overwrite,ref_pdb_fname,mark_resis,highlight,highlight_mode,skip_resis,ticks,dpi,title=None):
    if title is None:
        title = struct+' '+times+' resi '+plot_type+' correl ('+detail_level+', '+cmap+')'
        # New HAC-approved title format
        # TODO: FIXME: BAD HACK
        long_struct_names = {'1rx1':'closed-loop',
                             '1ra1':'open-loop',}
        _struct=long_struct_names[struct]
        _times = float(times[2:])
        _times = _times - 1.5 # equilibration
        _times = '%05.2f ns - %05.2f ns'%(times - 0.5, times + 0.5)
        title = _struct + ', ' + _times + ', ' + plot_type + ' correlation'
        #TODO: FIXME: end of BAD HACK
    figname = os.path.join(out_dir,title+'.png')
    if os.path.exists(figname) and not overwrite:
        sys.stdout.write("Skipping "+figname+'\n')
    step = {'fine':0.01,'coarse':0.05}[detail_level]
    #
    # Now we split up the mark_resis. If it's just a list of resis,
    # we deal with it as such. If it's two lists, we'll use the
    # first list on the left axis and the second on the bottom one.
    #
    if mark_resis and type(mark_resis[0]) in (type(()),type([])) and len(mark_resis) == 2:
        mark_left_resis,mark_bottom_resis = mark_resis
        mark_resis = sorted(mark_left_resis + mark_bottom_resis)
    else:
        mark_left_resis,mark_bottom_resis = mark_resis,mark_resis
    if plot_type in 'straight mainheavy allheavy sidechainhbond hbond'.split():
        #
        # These are the per-atom ones
        #
        master_residue_list,reference_coords,ca_atom_id_list,hydro_atom_id_list,nonhydro_atom_id_list,nosp_atom_id_list,nonnosp_atom_id_list,mainchain_atom_id_list,sidechain_atom_id_list,mainchain_nonhydro_atom_id_list,sidechain_hydro_atom_id_list,sidechain_nonhydro_atom_id_list,skip_atom_id_list = md_analysis_utils.setup_ca(ref_pdb_fname,indexing=0,skip_resis=skip_resis)
        # We cannot to mark the before we yank stuff out.
        # If we do, we change the size of the array, and end up# taking the wrong
        # stuff out.  So, we set things up here and mark later.
        # That's also why left_bar must be padded later on.
        atom_id_map = md_analysis_utils.get_resi_to_atom_id_map(ref_pdb_fname,indexing=0)
        mark_left_atoms = []
        mark_bottom_atoms = []
        for resi in mark_left_resis:
            mark_left_atoms += atom_id_map[resi]
        for resi in mark_bottom_resis:
            mark_bottom_atoms += atom_id_map[resi]
        orig_size = data.shape[0]
        bottom_bar = N.array([[(sentinel1 if i in mark_bottom_atoms else sentinel2) for i in range(orig_size)]])
        left_bar =   N.array([(sentinel1 if i in mark_left_atoms else sentinel2) for i in range(orig_size)])
        left_bar.shape = (len(left_bar),1)
        
        if plot_type == 'mainheavy':
            data = scipy.take(data, mainchain_nonhydro_atom_id_list, axis = 0)
            left_bar = scipy.take(left_bar, mainchain_nonhydro_atom_id_list, axis = 0)
            data = scipy.take(data, mainchain_nonhydro_atom_id_list, axis = 1)
            bottom_bar = scipy.take(bottom_bar,mainchain_nonhydro_atom_id_list, axis = 1)
        elif plot_type == 'allheavy':
                totake = mainchain_nonhydro_atom_id_list+sidechain_nonhydro_atom_id_list
                totake.sort()
                data = scipy.take(data, totake, axis = 0)
                left_bar = scipy.take(left_bar, totake, axis = 0)
                data = scipy.take(data, totake, axis = 1)
                bottom_bar = scipy.take(bottom_bar,totake, axis = 1)
            
                
        elif plot_type == 'sidechainhbond':
            data = scipy.take(data,sidechain_hydro_atom_id_list, axis = 0)
            left_bar = scipy.take(left_bar,sidechain_hydro_atom_id_list, axis = 0)
            data = scipy.take(data,[i for i in nosp_atom_id_list if i not in mainchain_nonhydro_atom_id_list], axis = 1)
            bottom_bar = scipy.take(bottom_bar,[i for i in nosp_atom_id_list if i not in mainchain_nonhydro_atom_id_list], axis = 1)
            if 0 in data.shape:
                print "COULD NOT CALCULATE sidechainhbond"
                return
            #data = scipy.take(data,nosp_atom_id_list,    axis = 1)
        elif plot_type == 'hbond':
            # We could also consider a masked array.
            # We could use for i in atom_id_list:data[...,i] = sentinel
            # or a real mask.
            data       = scipy.take(data,       hydro_atom_id_list, axis = 0)
            left_bar   = scipy.take(left_bar,   hydro_atom_id_list, axis = 0)
            data       = scipy.take(data,       nosp_atom_id_list , axis = 1)
            bottom_bar = scipy.take(bottom_bar, nosp_atom_id_list , axis = 1)
            if 0 in data.shape:
                print "COULD NOT CALCULATE hbond"
                return
    else:
        #
        # These are the per-residue ones
        #
        data = scipy.take(data,[i for i in range(data.shape[0]) if i not in skip_resis],axis=0)
        data = scipy.take(data,[i for i in range(data.shape[1]) if i not in skip_resis],axis=1)
        orig_size = data.shape[0]
        bottom_bar = N.array([[(sentinel1 if i in mark_bottom_resis else sentinel2) for i in range(orig_size)]])
        left_bar = [(sentinel1 if i in mark_left_resis else sentinel2) for i in range(orig_size)]
        new_shape = (len(left_bar),1)
        left_bar = N.array(left_bar)
        left_bar.shape = new_shape

    if mark_resis:
        #
        # When we mark the matrix, we make the bottom bar the size of the original
        # data and we pad the left_bar with enough slack to take care of the bottom-
        # left corner.  After we're done marking, we will use the bottom bar and the
        # left bar for highlighting.  So, we'll need to go ahead and pad the bottom
        # bar at that point.
        #
        num_blocks = int(0.01*max(data.shape))+1
        for i in range(num_blocks): data = N.concatenate((bottom_bar,data))
        padding = N.array([sentinel2,]*num_blocks)
        padding.shape = (num_blocks,1)
        new_shape = (len(left_bar)+len(padding),1)
        left_bar = N.concatenate((padding,left_bar))
        for i in range(num_blocks): data = N.concatenate((left_bar,data),axis=1)
        padding.shape = 1,num_blocks
        bottom_bar = N.concatenate((padding,bottom_bar),axis=1)

            
            
        
    _cmap,_norm = sentinel_maps_and_norms[cmap]
    #plotter = pylab.pcolormesh
    plotter = pylab.imshow
    print "I will plot now!"
    plotterargs = {pylab.imshow:{'interpolation':'nearest','origin':'lower'},
                   pylab.pcolormesh:{'shading':'flat'}}
    
    a = plotter(data,cmap=_cmap,norm=_norm,
                vmin=-1.0,vmax=1.0,**plotterargs[plotter])
    # For some reason, it draws the axes a little too long.
    if plot_type not in 'sidechainhbond hbond'.split():
        pylab.axis((0,data.shape[0],0,data.shape[1]))
    pylab.colorbar()
    pylab.title(title)
    if mark_resis and highlight:
        # Highlight the selected parts
        if highlight_mode in ('positive','negative'):
            mask = N.zeros(data.shape)
            for i in range(data.shape[0]):
                if left_bar[i] == sentinel1:
                    mask[i,] = True
            for j in range(data.shape[1]):
                if bottom_bar[0,j] == sentinel1:
                    mask[:,j] = True
            if highlight_mode == 'positive':
                data = N.ma.masked_where(N.logical_not(mask),data)
            elif highlight_mode == 'negative':
                data = N.ma.masked_where(mask,data)
            _cmap.set_bad('white',alpha=highlight)
            plotter(data,cmap=_cmap,norm=_norm,
                    vmin=-1.0,vmax=1.0,**plotterargs[plotter])
            
        elif highlight_mode == 'supernegative':
            last_val = left_bar[0]
            changes = []
            for i,val in enumerate(left_bar):
                if i == 0: continue
                if val != last_val:
                    changes.append(i)
                last_val = val
            for start,stop in zip(changes[:-1:2],changes[1::2]):
                facecolor='white'
                edgecolor='white'
                linewidth=0.0
                pylab.axvspan(start,stop,0.01,1,fc=facecolor,ec=edgecolor,lw=linewidth,alpha=highlight)
                pylab.axhspan(start,stop,0.01,1,fc=facecolor,ec=edgecolor,lw=linewidth,alpha=highlight)
            for i in changes:
                # I don't actually like the way this looks
                continue
                pylab.axhline(y=i,color='black')
                pylab.axvline(x=i,color='black')
    if ticks:
        from matplotlib.ticker import MultipleLocator, FormatStrFormatter, Formatter, FixedLocator, LinearLocator, Locator

        class OffsetFormatStrFormatter(Formatter):
            def __init__(self,fmt,offset):
                self.fmt = fmt
                self.offset = offset
            def __call__(self,x,pos=None):
                return self.fmt % (x-self.offset)

        class ExplicitLinearLocator(Locator):
            def __init__(self,vmin,vmax,numticks):
                self.vmin = vmin
                self.vmax = vmax
                self.numticks=numticks
            def __call__(self):
                'Return the location of the ticks'
                ticklocs = mlab.linspace(self.vmin,self.vmax,self.numticks)
                return ticklocs

        offset = (num_blocks if mark_resis else 0)
        ymajorFormatter = OffsetFormatStrFormatter('%d',offset)
        xmajorFormatter = OffsetFormatStrFormatter('%d',offset)
        yminorLocator   = ExplicitLinearLocator(offset,data.shape[0],81)
        xminorLocator   = ExplicitLinearLocator(offset,data.shape[1],81)

        ax = pylab.gca()
        ax.yaxis.set_major_formatter(ymajorFormatter)
        ax.yaxis.set_minor_locator(yminorLocator)
        pylab.yticks(mlab.linspace(offset,data.shape[0],9))

        ax.xaxis.set_major_formatter(xmajorFormatter)
        ax.xaxis.set_minor_locator(xminorLocator)
        pylab.xticks(mlab.linspace(offset,data.shape[1],9))
    else:
        pylab.xticks([])
        pylab.yticks([])
                 
    if save_fig:
        sys.stdout.write("Saving "+title+'\n')
        sys.stdout.flush()
        pylab.savefig(figname,dpi=dpi)
        pylab.clf()
    else:
        pylab.show()
        print "now showing!"
        if 1: pylab.savefig('./x',dpi=dpi)
    del data

def make_correl_plots_for_movie(structures=['1ra1','1rx1'],
                                all_times=['NS1.0,'],
                                cmaps=['Normal',],
                                detail_levels=['fine',],
                                plot_types='abs max min avg straight mainheavy allheavy sidechainhbond hbond ca'.split(),
                                overwrite=True,
                                image_dir='BigMovieImages/',
                                dat_dir='/Users/mglerner/work/Dynamics-DHFR/BigDynamicsMoviePtrajFiles/',
                                ref_pdb_fname = '/Users/mglerner/work/Dynamics-DHFR/MD_Files/1RX1/ChainAPDBFilesFromTrajectory/AlignedVs1RX1/1rx1_trajectory_A_snap_0001_trans.pdb',
                                mark_resis=[],
                                highlight=0.0,
                                highlight_mode='positive',
                                skip_resis=[],
                                ticks=True,
                                dpi=200,
                                title=None,
                                ):
    '''
    This will generate all of the .png images that you need to make your
    movies.  It will also spit out an html file called
    BigMovieImages<something>.html for you to look at.
    
    Empirical data: Setting width="86%" on all of the <img/> tags will make it so that you can
    nicely fit the title and four rows of pictures on one page.
    '''
    
    for cmap in cmaps:
        for detail_level in detail_levels:
            fname = 'BigMovieImages'+cmap+detail_level+'.html'
            fout = file(fname,'w')
            fout.write('<html><head><title>Big Movie Images</title></head><body>')
            fout.write('<h1>'+cmap+' '+detail_level+'</h1><table>')
            made_plots_already = False
            for plot_type in plot_types:
                for times in all_times:
                    fout.write('<tr>')
                    for struct in structures:
                        title = struct+' '+times+' resi '+plot_type+' correl ('+detail_level+', '+cmap+')'
                        # New HAC-approved title format
                        # TODO: FIXME: BAD HACK
                        long_struct_names = {'1rx1':'closed-loop',
                                             '1ra1':'open-loop',}
                        _struct=long_struct_names[struct]
                        _times = float(times[2:])
                        _times = _times - 1.5 # equilibration
                        _times = '%05.2f ns - %05.2f ns'%(_times - 0.5, _times + 0.5)
                        title = _struct + ', ' + _times + ', ' + plot_type + ' correlation'
                        #TODO: FIXME: end of BAD HACK
                        image_name = os.path.join(image_dir,title+'.png')
                        if os.path.exists(image_name) and not made_plots_already:
                            if overwrite:
                                print "Overwriting",title
                        if not made_plots_already:
                            make_several_correl_plots(struct,times,cmap,detail_level,save_fig=True,out_dir=image_dir,overwrite=overwrite,
                                                      dat_dir=dat_dir,
                                                      ref_pdb_fname = ref_pdb_fname,
                                                      plot_types=plot_types,
                                                      mark_resis=mark_resis,
                                                      highlight=highlight,
                                                      highlight_mode=highlight_mode,
                                                      skip_resis=skip_resis,
                                                      ticks=ticks,
                                                      dpi=dpi,
                                                      title=title,
                                                      )
                            
                        fout.write('<td><img src="'+image_name+'" width="86%" /></td>')
                    fout.write('</tr>')
                    fout.flush()
                fout.write('</tr>')
                fout.flush()
                
                made_plots_already=True
            fout.write('</table>')
            fout.write('</body></html>')
            fout.close()
    
def plot_rmsds_for_DHFR_paper():
    #start_time = 0.505 # in NS
    start_time = 0 # in NS
    stop_time  = 10.505 # in NS
    for time in 'xtal 505ps 1505ps 3005ps'.split():
        pylab.clf()
        data = {}
        data_dir = '/Users/mglerner/work/Dynamics-DHFR/MD_Files/RMSDAnalysis'
        data['rx_all'] = scipy.io.read_array(os.path.join(data_dir,'1rx1_%s_all_atom_rmsds.txt'%time))
        data['rx_ca']  = scipy.io.read_array(os.path.join(data_dir,'1rx1_%s_ca_rmsds.txt'%time))
        data['ra_all'] = scipy.io.read_array(os.path.join(data_dir,'1ra1_%s_all_atom_rmsds.txt'%time))
        data['ra_ca']  = scipy.io.read_array(os.path.join(data_dir,'1ra1_%s_ca_rmsds.txt'%time))
        #for struct in 'ra_all rx_all ra_ca rx_ca'.split():
        for struct in 'ra_all rx_all'.split():
        #for struct in 'ra_ca rx_ca'.split():
            labels = {#'ra_all':'1RA1 heavy atom RMSD vs %s (Avg. after = %s)',
                      'ra_all':'Open loop (Avg. = %s)',
                      'ra_ca' :'1RA1 ca RMSD vs %s                (Avg. after = %s)',
                      #'rx_all':'1RX1 heavy atom RMSD vs %s (Avg. after = %s)',
                      'rx_all':'Closed loop (Avg. = %s)',
                      'rx_ca' :'1RX1 ca RMSD vs %s                (Avg. after = %s)',
                      }
            this_data = [(i[0]*5/1000.0,i[1]) for i in data[struct] if ((i[0]*5/1000.0) <= stop_time) and ((i[0]*5/1000.0) >= start_time)]
            data_after = [i[1] for i in this_data if (i[0] >= {'xtal':0,'505ps':0.505,'1505ps':1.505,'3005ps':3.005}[time])]
            avg_after = sum(data_after)/len(data_after)
            print struct,avg_after,len(data_after)
            #pylab.plot([i[0]for i in this_data],[i[1] for i in this_data],label=labels[struct]%(time,avg_after))
            pylab.plot([i[0]for i in this_data],[i[1] for i in this_data],label=labels[struct]%(avg_after))
            #pylab.plot([i[0]for i in this_data],[avg_after for i in this_data],color="black")
            f = file('/Users/mglerner/tmp/%s_MD_RMSDs_vs_%s.csv'%(struct,time),'w')
            f.write('Time,RMSD,Average %s RMSD after %s = %s\n'%(struct,time,avg_after))
            for i in this_data:
                f.write('%s,%s\n'%i)
            f.close()
        close_ticks = False
        if close_ticks:
            pylab.xticks([i/2.0 for i in range(24)])
            pylab.yticks([i/10.0 for i in range(30)])
            pylab.axis((start_time,stop_time,0,3))
            pylab.legend()
            pylab.grid(b=True)
            pylab.savefig(os.path.join(data_dir,'1RA1and1RX1_MD_RMSDs_vs_%s_just_ca.png'%time))
        else:
            pylab.xticks([i/2.0 for i in range(24)])
            pylab.yticks([i/2.0 for i in range(30)])
            pylab.axis((start_time,stop_time,0,3))
            pylab.legend()
            pylab.grid(b=True)
            pylab.title('Closed- and Open-loop RMSDs vs. %s structure'%{'xtal':'crystal'}.get(time,time))
            a = input()
            pylab.savefig(os.path.join(data_dir,'1RA1and1RX1_MD_RMSDs_vs_%s_just_ca_fewer_ticks.png'%time))
            
    
def plot_Sander_Data_for_DHFR_paper():
    pylab.clf()
    all_types = '1_4_EEL EELEC TEMP_K_ 1_4_NB EHBOND VDWAALS ANGLE EKCMT VIRIAL BOND EKtot VOLUME CONSTRAINT EPtot DIHED Etot Density NSTEP EAMBER__non_constraint_ PRESS'.split()
    pos_energy_types = 'EKtot VDWAALS VIRIAL EKCMT 1_4_EEL ANGLE DIHED 1_4_NB BOND'.split()
    neg_energy_types = 'Etot EPtot EELEC'.split()
    energy_types = pos_energy_types + neg_energy_types
    non_energies = 'VOLUME Density PRESS TEMP_K_'.split()
    for struct in '1rx1'.split():
    #for struct in '1ra1 1rx1'.split():
        start_dir = '/Users/mglerner/work/Dynamics-DHFR/MD_Files/AllSanderOutFiles/Combined/%sfiles/'%struct
        if 0:
            # one-at-a-time
            for et in energy_types:
                data = scipy.io.read_array(os.path.join(start_dir,et+'.txt'))
                pylab.plot([(i[0]-500)/1000.0 for i in data],[i[1] for i in data],color='gray',label=struct +'    '+et.replace('_',' '))
                last_part_of_data = [i[1] for i in data if ((i[0]-500)/1000.0) >= 1.5]
                avg = sum(last_part_of_data)/float(len(last_part_of_data))
                pylab.plot([(i[0]-500)/1000.0 for i in data],[avg for i in data],color='black',label='Average after 1.5ns (%.3g)'%avg)
                print et,'\t %.3g'%avg
                pylab.xticks([i/2.0 for i in range(24)])
                pylab.grid(b=True)
                a = pylab.axis()
                pylab.axis((0,5.5,a[2],a[3]))
                pylab.legend(loc=0)
                pylab.savefig(os.path.join('/Users/mglerner/work/DHFR/MD_Analysis',struct+'_'+et+'.png'))
                pylab.clf()
        if 1:
            # all at once
            for et in pos_energy_types:
                data = scipy.io.read_array(os.path.join(start_dir,et+'.txt'))
                avg = sum([i[1] for i in data])/float(len(data))
                pylab.plot([(i[0]-500)/1000.0 for i in data],[i[1] for i in data],label=et+'(%.3g)'%avg)
                print et,'\t %.3g'%avg
                pylab.legend()
            pylab.xticks([i/2.0 for i in range(24)])
            pylab.grid(b=True)
            pylab.legend(loc=0)
            pylab.show()
        
    
def make_animated_gif():
    thumbPart = '-size 200x200 -geometry 200x200'
    for k in fileMap.keys():
        ps = getFilename(k,'.ps',imgDir)
        gif = getFilename(k,'.gif',imgDir)
        thumb = getFilename(k,'_thumb.gif',imgDir)
        os.system('/usr/bin/env convert -rotate 90 %s %s' % (ps,gif))
        os.system('/usr/bin/env convert -rotate 90 %s %s %s' % (thumbPart,
                                                                ps,
                                                                thumb))
    
        


if __name__ == '__main__':
    #plotrms(glob.glob('rms*.txt'))
    #plotrms2()
    #plot_ile_leu_dist()
    #plotcontours()
    
    #plotrms3()
    #plotcontours_in_several_ways()
    #plot_aligned_rms()
    #plot_ileleu()

    #plot_diff_colormesh('/Users/mglerner/work/Dynamics-DHFR/MD_Files/1RX1/1rx1_byres_correlmat.dat','/Users/mglerner/work/Dynamics-DHFR/MD_Files/1RA1/1ra1_byres_correlmat.dat')

    if 0:
        plot_rmsds_for_DHFR_paper()
    elif 0:
        #
        # Current correl/covar plotting mechanisms
        #
        #make_all_correl_plots()
        make_correl_plots_for_movie()
    elif 0:
        #
        # Current RMS plotting mechanisms for pymol-generated
        #
        pymol_generated_structure_parts = 'all_atom_substrate_binding_full all_atom_substrate_binding_min n_substrate_binding_full n_substrate_binding_min ca_substrate_binding_full ca_substrate_binding_min all_atom_substrate_binding_min_noR52 n_substrate_binding_min_noR52 ca_substrate_binding_min_noR52'.split()
        for (out_dir,structure_parts) in (('Standard',[i for i in pymol_generated_structure_parts if 'R52' in i]),
                                          ):
            for structure_part in structure_parts:
                for struct in '1rx1 1ra1'.split():
                    for include_averages in (True,):#,False):
                        plot_aligned_rms(struct=struct,
                                         structure_part=structure_part,
                                         out_dir = out_dir,
                                         include_averages=include_averages,
                                         prefix='pymol_rmsd',
                                         )
    elif 0:
        #
        # Current RMS plotting mechanisms
        #
        #for structure_part in 'm20 FG CD GH all'.split():
        standard_parts = 'm20 all FG CD GH noloops subdomain1 subdomain2 substrate_binding substrate_binding2'.split() 
        residue_chunks = '0-20 20-40 40-60 60-80 80-100 100-120 120-140 140-160'.split()
        anti_correlated_parts = 'm20 40-50 59-68 71-74 95-97'.split() + '116-124 142-149'.split() # 116-124 is correlated with 40-50 and 59-68, 142-149 is correlated with 40-50
        correlated_parts = '110-114 134-141 107-113 149-156 132-142 144-157'.split() + '115-125 1-10 58-62 43-49 38-63 42-51'.split()
        for (out_dir,structure_parts) in (#('Standard',standard_parts),
                                          #('ResidueChunks',residue_chunks),
                                          #('AntiCorrelated',anti_correlated_parts),
                                          #('Correlated',correlated_parts),
                                          ('Standard',[i for i in standard_parts if i.startswith('substrate_binding')]),
                                          ):
            for structure_part in structure_parts:
                for struct in '1ra1 1rx1'.split():
                    for include_averages in (True,):#,False):
                        plot_aligned_rms(struct=struct,
                                         structure_part=structure_part,
                                         out_dir = out_dir,
                                         include_averages=include_averages,
                                         prefix='rmsd_ca',
                                         )
        
