#!/usr/bin/env python
import numpy as npy
from os import listdir
import os.path
from re import match
from time import time
import cPickle

def savespkrate(prefix,dn):
    if len(prefix)>0: 
        outname = prefix+'all.pkl'
    fns = [fn for fn in listdir(dn) if fn.endswith('.pkl') and fn.startswith(prefix) and 'run' in fn]
    f = []
    for fn in fns:
        f.append(cPickle.load(open(os.path.join(dn,fn))))
    f = npy.array(f)
    cPickle.dump(f,open(os.path.join(dn,outname),'w'))

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
        m = match('(.*)run\d+.pkl',l)
        if m and not os.path.exists(os.path.join(dn,'%sall.pkl' % (m.groups()[0]))):
            prefixes.append(m.groups()[0])

    prefixes = list(set(prefixes))
    for runnum, prefix in enumerate(prefixes): 
        runsavespkrate(prefix,dn)

updatedir('data')
