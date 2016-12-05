#!/usr/bin/env python
"""

TODO:

1) Get rid of unused functions
2) Possibly move the cmap lookup to tool_utils.

"""

from .sentinel_map import SentinelMap,SentinelNorm
import glob
#import matplotlib.numerix.ma as ma
from . import md_analysis_utils
import sys
import pylab,matplotlib,os
from matplotlib import mlab
import numpy as np
import scipy.io
import pprint
import bz2
import copy
from .tool_utils import read_data
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
sentinel_maps_and_norms = {'Normal':(SentinelMap(my_cmap),#,sentinels={sentinel1:rgb1,sentinel2:rgb2}),
                                     SentinelNorm(vmin=-1.0,vmax=1.0)#,ignore=[sentinel1,sentinel2,])
                                     ),
                           'Scaled':(SentinelMap(my_scaled_cmap),#sentinels={sentinel1:rgb1,sentinel2:rgb2}),
                                     SentinelNorm(vmin=-1.0,vmax=1.0)#,ignore=[sentinel1,sentinel2,])
                                     ),
                           }
    
                                    
##                           'jet'   :(SentinelMap(pylab.cm.jet,sentinels={sentinel1:rgb1,sentinel2:rgb2}),
##                                     SentinelNorm(ignore=[sentinel1,sentinel2,],vmin=-1.0,vmax=1.0)
##                                     ),
##                           }
for map in [m for m in list(pylab.cm.datad.keys()) if not m.endswith("_r")]:
    sentinel_maps_and_norms[map] = (SentinelMap(pylab.cm.get_cmap(map)),#sentinels={sentinel1:rgb1,sentinel2:rgb2}),
                                     SentinelNorm(vmin=-1.0,vmax=1.0)#,ignore=[sentinel1,sentinel2,])
                                     #SentinelNorm(ignore=[sentinel1,sentinel2,],vmin=-1.0,vmax=1.0)
                                    )

                                     

def make_several_correl_plots(struct,
                              times,
                              cmap,
                              detail_level,
                              save_fig,
                              out_dir,
                              overwrite,
                              dat_dir,
                              ref_pdb_fname,
                              plot_types,#='ca avg max min abs straight mainheavy allheavy sidechainhbond hbond'.split(),
                              mark_resis,
                              highlight,
                              highlight_mode,
                              skip_resis,
                              ticks,
                              dpi,
                              titleprefix,
                              ):
    for plot_type in [i for i in 'abs max min avg ca'.split() if i in plot_types]:
        if plot_type == 'avg': fname = os.path.join(dat_dir,'%s/%s_%s_byres_correlmat.dat'%(struct,struct,times))
        else:                  fname = os.path.join(dat_dir,'%s/%s_%s_all_atom_correlmat_resi_%s.dat'%(struct,struct,times,plot_type))
        data = read_data(fname)
        if data is None:
            continue
        title = titleprefix +plot_type + ' correlation'
        print(title,'S')
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
            title = titleprefix +plot_type + ' correlation'
            print(title,'L')
            _make_one_correl_plot(data,struct,times,cmap,detail_level,plot_type,save_fig,dat_dir,out_dir,overwrite,ref_pdb_fname,mark_resis,highlight,highlight_mode,skip_resis,ticks,dpi,title)
            del _data
        del data


def _make_one_correl_plot(data,struct,times,cmap,detail_level,plot_type,save_fig,dat_dir,out_dir,overwrite,ref_pdb_fname,mark_resis,highlight,highlight_mode,skip_resis,ticks,dpi,title):

    print("Title is",title)
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
        bottom_bar = np.array([[(sentinel1 if i in mark_bottom_atoms else sentinel2) for i in range(orig_size)]])
        left_bar =   np.array([(sentinel1 if i in mark_left_atoms else sentinel2) for i in range(orig_size)])
        left_bar.shape = (len(left_bar),1)
        
        if plot_type == 'mainheavy':
            data = scipy.take(data, mainchain_nonhydro_atom_id_list, axis = 0)
            left_bar = scipy.take(left_bar, mainchain_nonhydro_atom_id_list, axis = 0)
            data = scipy.take(data, mainchain_nonhydro_atom_id_list, axis = 1)
            bottom_bar = scipy.take(bottom_bar,mainchain_nonhydro_atom_id_list, axis = 1)
        elif plot_type == 'allheavy':
                totake = sorted(mainchain_nonhydro_atom_id_list+sidechain_nonhydro_atom_id_list)
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
                print("COULD NOT CALCULATE sidechainhbond")
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
                print("COULD NOT CALCULATE hbond")
                return
    else:
        #
        # These are the per-residue ones
        #
        skip_resis = [i-1 for i in skip_resis]
        data = scipy.take(data,[i for i in range(data.shape[0]) if i not in skip_resis],axis=0)
        data = scipy.take(data,[i for i in range(data.shape[1]) if i not in skip_resis],axis=1)
        orig_size = data.shape[0]
        bottom_bar = np.array([[(sentinel1 if i in mark_bottom_resis else sentinel2) for i in range(orig_size)]])
        left_bar = [(sentinel1 if i in mark_left_resis else sentinel2) for i in range(orig_size)]
        new_shape = (len(left_bar),1)
        left_bar = np.array(left_bar)
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
        for i in range(num_blocks): data = np.concatenate((bottom_bar,data))
        padding = np.array([sentinel2,]*num_blocks)
        padding.shape = (num_blocks,1)
        new_shape = (len(left_bar)+len(padding),1)
        left_bar = np.concatenate((padding,left_bar))
        for i in range(num_blocks): data = np.concatenate((left_bar,data),axis=1)
        padding.shape = 1,num_blocks
        bottom_bar = np.concatenate((padding,bottom_bar),axis=1)

       
    _cmap,_norm = sentinel_maps_and_norms[cmap]
    #plotter = pylab.pcolormesh
    plotter = pylab.imshow
    print("I will plot now!")
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
            mask = np.zeros(data.shape)
            for i in range(data.shape[0]):
                if left_bar[i] == sentinel1:
                    mask[i,] = True
            for j in range(data.shape[1]):
                if bottom_bar[0,j] == sentinel1:
                    mask[:,j] = True
            if highlight_mode == 'positive':
                data = np.ma.masked_where(np.logical_not(mask),data)
            elif highlight_mode == 'negative':
                data = np.ma.masked_where(mask,data)
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
                ticklocs = np.linspace(self.vmin,self.vmax,self.numticks)
                return ticklocs

        offset = (num_blocks if mark_resis else 0)
        ymajorFormatter = OffsetFormatStrFormatter('%d',offset)
        xmajorFormatter = OffsetFormatStrFormatter('%d',offset)
        yminorLocator   = ExplicitLinearLocator(offset,data.shape[0],81)
        xminorLocator   = ExplicitLinearLocator(offset,data.shape[1],81)

        ax = pylab.gca()
        ax.yaxis.set_major_formatter(ymajorFormatter)
        ax.yaxis.set_minor_locator(yminorLocator)
        pylab.yticks(np.linspace(offset,data.shape[0],9))

        ax.xaxis.set_major_formatter(xmajorFormatter)
        ax.xaxis.set_minor_locator(xminorLocator)
        pylab.xticks(np.linspace(offset,data.shape[1],9))
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
        print("now showing!")
        if 1: pylab.savefig('./x',dpi=dpi)
    del data

def make_correl_plots_for_movie(structures,
                                all_times,
                                window_size_ns,
                                cmaps,
                                detail_levels,
                                plot_types,#='abs max min avg straight mainheavy allheavy sidechainhbond hbond ca'.split(),
                                overwrite,
                                image_dir,
                                dat_dir,
                                ref_pdb_fname,
                                mark_resis,
                                highlight,
                                highlight_mode,
                                skip_resis,
                                ticks,
                                dpi,
                                title,
                                ):
    '''
    This will generate all of the .png images that you need to make your
    movies.  It will also spit out an html file called
    BigMovieImages<something>.html for you to look at.
    
    Empirical data: Setting width="86%" on all of the <img/> tags will make it so that you can
    nicely fit the title and four rows of pictures on one page.

    plot_types can be 'abs max min avg straight mainheavy allheavy sidechainhbond hbond ca'.split(),
    '''
    for cmap in cmaps:
        for detail_level in detail_levels:
            fname = 'BigMovieImages'+cmap+detail_level+'.html'
            fout = open(fname,'w')
            fout.write('<html><head><title>Big Movie Images</title></head><body>')
            fout.write('<h1>'+cmap+' '+detail_level+'</h1><table>')
            made_plots_already = False
            for plot_type in plot_types:
                # NOTE: This loop only happens once. Then we set
                #       "made_plots_already" to true We actually make
                #       all of the plots in the first call to
                #       make_several_correl_plots. This trickery is
                #       just so that it's really easy to write out the
                #       HTML file.
                #
                #       For the most part, this is easy to deal with.
                #       One subtlely is with the titles. We calculate
                #       them here to stamp them into the HTML files,
                #       but we must also calculate them in the
                #       make_several_correl_plots function.
                for times in all_times:
                    fout.write('<tr>')
                    for struct in structures:
                        title = struct+' '+times+' resi '+plot_type+' correl ('+detail_level+', '+cmap+')'
                        # New HAC-approved title format
                        # TODO: FIXME: BAD HACK
                        long_struct_names = {'1rx1':'closed-loop',
                                             '1ra1':'open-loop',}
                        _struct=long_struct_names.get(struct,struct)
                        _times = float(times[2:])
                        _times = '%05.2f ns - %05.2f ns'%(_times - window_size_ns/2., _times + window_size_ns/2.)
                        titleprefix = _struct + ', ' + _times + ', ' 
                        title = titleprefix + plot_type + ' correlation'
                        print(title)
                        image_name = os.path.join(image_dir,title+'.png')
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
                                                      titleprefix=titleprefix,
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
    
    
        


