#!/usr/bin/env python
import numpy as npy
import tables as pyt
from os import listdir
import os.path
from re import match
from time import time
from common import h5get

filters = pyt.Filters(complevel=9, complib='zlib', shuffle=True, fletcher32=True)

dt = .1/1000  # in seconds
spkth = -20

def savespkrate(prefix,dn,tstop=None,dt=dt,outname='some-spkavg.h5'):
    if len(prefix)>0: 
        outname = prefix+'spkavg.h5'
    h5fns = [fn for fn in listdir(dn) if fn.endswith('.h5') and fn.startswith(prefix) and 'run' in fn]
    shape = h5get(os.path.join(dn,h5fns[0]),'s.shape')
    matintraspk = npy.zeros(shape[:-1]+(len(h5fns),))
    matintraerr = npy.zeros(shape[:-1]+(len(h5fns),))
    matintramus = npy.zeros(shape[:-1]+(len(h5fns),))
    for ind,h5fn in enumerate(h5fns):
        shape = list(h5get(os.path.join(dn,h5fn),'s.shape'))
        tslice = npy.arange(shape[-1])
        shape[-1] = tslice.size
        tspan = shape[-1]*dt
        trace = h5get(os.path.join(dn,h5fn),'s')
        matintraspk[...,ind] = (npy.diff(npy.array(trace[...,tslice]>spkth,npy.int32),axis=-1)==1).sum(axis=-1)/tspan
        matintraerr[...,ind] = npy.std(trace[...,tslice],axis=-1)
        matintramus[...,ind] = npy.mean(trace[...,tslice],axis=-1)

    matinterspk = matintraspk.mean(axis=len(shape)-1)
    matintererr = matintraspk.std(axis=len(shape)-1)
    matintermus = matintramus.mean(axis=len(shape)-1)
    h5f = pyt.openFile(os.path.join(dn,outname),'w',filters=filters)
    h5f.createCArray(h5f.root, 'spikes', shape=matinterspk.shape, atom=pyt.Float64Atom(), chunkshape=matinterspk.shape)[:] = matinterspk
    h5f.createCArray(h5f.root, 'spikesstd', shape=matintererr.shape, atom=pyt.Float64Atom(), chunkshape=matintererr.shape)[:] = matintererr
    h5f.createCArray(h5f.root, 'spikesindv', shape=matintraspk.shape, atom=pyt.Float64Atom(), chunkshape=matintraspk.shape)[:] = matintraspk
    h5f.createCArray(h5f.root, 'spikesindvstd', shape=matintraerr.shape, atom=pyt.Float64Atom(), chunkshape=matintraerr.shape)[:] = matintraerr
    h5f.createCArray(h5f.root, 'means', shape=matintermus.shape, atom=pyt.Float64Atom(), chunkshape=matintermus.shape)[:] = matintermus
    h5f.createCArray(h5f.root, 'meansindv', shape=matintramus.shape, atom=pyt.Float64Atom(), chunkshape=matintramus.shape)[:] = matintramus
    h5f.close()

def runsavespkrate(prefix,dn):
    print 'prefix: %s' % prefix
    try:
        tic = time()
        savespkrate(prefix,dn)
        print '%s took %d seconds to process' % (prefix,time()-tic)
    except:
        print "%s did not get processed!" % prefix
        pass

def updatedir(dn='.'):
    ls = sorted(listdir(dn))
    prefixes = []
    for l in ls:
        m = match('(.*)run\d+.h5',l)
        if m and not os.path.exists(os.path.join(dn,'%sspkavg.h5' % (m.groups()[0]))): 
            prefixes.append(m.groups()[0])

    prefixes = list(set(prefixes))
    for runnum, prefix in enumerate(prefixes): 
        runsavespkrate(prefix,dn)

updatedir('data')
