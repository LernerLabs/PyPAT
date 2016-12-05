#!/usr/bin/env python2.4
#!/usr/bin/env python2.4
from __future__ import division
#from pylab import *
import pylab
from matplotlib.colors import Colormap, LinearSegmentedColormap, normalize
from matplotlib import numerix as nx
import matplotlib.numerix.ma as ma
from matplotlib.numerix import array, arange, alltrue
from types import IntType, FloatType, ListType

class SentinelMap(Colormap):
    def __init__(self, cmap=pylab.cm.jet, sentinels={}):
        # boilerplate stuff
        self.N = cmap.N
        self.name = 'SentinelMap'
        self.cmap = cmap
        self.sentinels = sentinels

        for rgb in sentinels.values():
            if len(rgb)!=3:
                raise ValueError('sentinel color must be RGB')
        self.is_gray = cmap.is_gray
        self.set_bad = cmap.set_bad


    def __call__(self, scaledImageData, alpha=1):
        # assumes the data is already normalized (ignoring sentinels)
        # clip to be on the safe side
        rgbaValues = self.cmap(nx.clip(scaledImageData, 0.,1.))
        for sentinel,rgb in self.sentinels.items():
            r,g,b = rgb
            if (scaledImageData==sentinel).max():
                rgbaValues[...,0] =  nx.where(scaledImageData==sentinel, r, rgbaValues[...,0])
                rgbaValues[...,1] =  nx.where(scaledImageData==sentinel, g, rgbaValues[...,1])
                rgbaValues[...,2] =  nx.where(scaledImageData==sentinel, b, rgbaValues[...,2])
                rgbaValues[...,3] =  nx.where(scaledImageData==sentinel, alpha, rgbaValues[...,3])
        return rgbaValues


class SentinelNorm(normalize):
    """
    Leave the sentinel unchanged
    """
    def __init__(self, ignore=[], vmin=None, vmax=None, clip = True):
        self.vmin=vmin
        self.vmax=vmax
        self.clip = clip
        
        if type(ignore) in [IntType, FloatType]:
            self.ignore = [ignore]
        else:
            self.ignore = list(ignore)
        self.ignore_mask=None
        
    def __call__(self, value, clip=None):
        
        if clip is None:
            clip = self.clip

        # ensure that we have a masked array val to work with
        if isinstance(value, (int, float)):
            vtype = 'scalar'
            val = ma.array([value])
        else:
            vtype = 'array'
            if ma.isMA(value):
                val = value
            else:
                val = ma.asarray(value)

        # create ignore_mask, val=sentinel1 | val= sentinel2..
        if self.ignore is not None:
            self.get_ignore_mask(val)
            
        # find min and max over points not masked by ignore_mask or by original mask of val
        self.autoscale(val)

        # now do scaling
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            if False in val.mask:
                raise ValueError("minvalue must be less than or equal to maxvalue")
            else:
                # array is completely masked. doesn't matter what values are for plot
                return 0.*value
        elif vmin==vmax:
            return 0.*value
        else:
            # scale points not masked by ignore_mask or by original mask of val
            scale = 1./(vmax-vmin)
            result = (val-vmin)*scale
            if clip:
                result = nx.clip(result,0.,1.)
            # set result over sentinel points to sentinel values 
            if self.ignore is not None:
                result[self.ignore_mask]=val.data[self.ignore_mask]
                
        if vtype == 'scalar':
            result = result[0]
        return result
    
    def get_ignore_mask(self, A):
        if ma.isMA(A):
            A=A.data
        if self.ignore is not None:
            self.ignore_mask = False
            for ignore in self.ignore:
                self.ignore_mask |= A==ignore


    def autoscale(self, A):
        # self.scaled is method in base class Normalize [colors.py], is True if self.vmin,vmax already defined
        if not self.scaled():
            if self.ignore is not None:
                if self.ignore_mask is None:
                    self.get_ignore_mask(A)
                A = ma.masked_where(self.ignore_mask,A)

            if self.vmin is None: self.vmin = A.min()
            if self.vmax is None: self.vmax = A.max()
            
    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = self.vmin, self.vmax

        if isinstance(value, (int, float)):
            return vmin + value * (vmax - vmin)
        else:
            val = ma.asarray(value)
            result = vmin + val * (vmax - vmin)
            if self.ignore is not None:
                if self.ignore_mask is None:
                    self.get_ignore_mask(value)
                result[self.ignore_mask]=val.data[self.ignore_mask]
            return result




cdict = {'red': ((0.0, 0.0, 0.0),
                 # This means that, at 0.5, Red starts at 1.0 and goes down
                 # to zero at 0.0.  It starts at 0.7 and goes up to 1.0 at 1.0.
                 (0.5, 1.0, 0.7), 
                 (1.0, 0.0, 0.0)
                 ),
         'green': (#(0.0, 0.0, 0.0),
                   #(0.5, 1.0, 0.0),
                   #(1.0, 1.0, 0.0)
                   (0.0,0.0,0.0),
                   (0.5,0.0,0.0),
                   (1.0,0.0,0.0)
                   ),
         'blue': (#(0.0, 0.0, 0.0),
                  #(0.5, 1.0, 0.0),
                  #(1.0, 0.5, 0.0)
                   (0.0,0.0,0.0),
                   (1.0,0.0,0.0)
                   )
         }
cdict = {'red':  (((-1.0+1)/2,32/255, 32/255),
                  ((-0.6+1)/2,32/255, 32/255),
                  ((-0.3+1)/2,32/255, 32/255),
                  ((-0.2+1)/2,0,0),
                  ((-0.1+1)/2,0,0),
                  (( 0.0+1)/2,0,0),
                  (( 0.1+1)/2,0,0),
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
if __name__ == '__main__':
    # define the sentinels
    sentinel1 = 10
    sentinel2 = -10
    # define the colormap and norm
    rgb1 = (0.,0.,0.)
    rgb2 = (1.,1.,1.)


    #my_cmap = SentinelLinearSegmentedColormap('my_colormap',cdict,256,sentinels={sentinel1:rgb1,sentinel2:rgb2})
    my_cmap = LinearSegmentedColormap('my_colormap',cdict,256)
    norm = SentinelNorm(ignore=[sentinel1,sentinel2,sentinel_zero,sentinel_zero2])
    my_sentinelcmap = SentinelMap(my_cmap,sentinels={sentinel1:rgb1,
                                                     sentinel2:rgb2,
                                                     })


    #pcolor(rand(10,10),cmap=my_cmap)
    n = 100
    X = array(pylab.outerproduct(arange(-1,1,0.01),ones(100)))
    # replace some data with sentinels
    X[int(.1*n):int(.2*n), int(.5*n):int(.7*n)]  = sentinel1
    X[int(.6*n):int(.9*n), int(.2*n):int(.5*n)]  = sentinel2

    # make mask
    mask = nx.mlab.rand(n*2,n) >0.5
    print 'mask\n ',mask
    print mask.shape
    print X.shape


    # now mask X
    X= ma.masked_where(mask,X)

    #pylab.pcolormesh(X,cmap=my_cmap,shading='flat')
    pylab.pcolormesh(X,cmap=my_sentinelcmap,shading='flat',norm=norm)
    pylab.colorbar()
    pylab.show()
