import matplotlib as mpl
mpl.use('Agg')  # for headless plotting
import tables as pyt
import numpy as npy
import pylab as pyl

# load data
def h5get(filename, arrayname, slice='[:]'):
    h5 = pyt.openFile(filename)
    h5dat = eval('h5.root.'+arrayname+slice)
    h5.close()
    return h5dat

# plotting functions
def setaspectsquare(ax, aspect=1.):
    ax.set_aspect(aspect*npy.ptp(ax.get_xlim())/npy.ptp(ax.get_ylim()))

def newfig(figsize=(8., 6.), dpi=91, num=None):
    h = pyl.figure(num, figsize=figsize, dpi=dpi)
    h.set_size_inches(figsize)
    return h
