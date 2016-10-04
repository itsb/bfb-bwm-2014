#!/usr/bin/env python
import sys
import cPickle
import numpy as npy

# numpy 'pretty print'
npy.set_printoptions(linewidth=135)
npy.set_printoptions(precision=3)
npy.set_printoptions(suppress=True)

from neuron import h
h.load_file('basal_project.hoc')

numbranches = int(sys.argv[1])
simiter = int(sys.argv[2])-1
if numbranches==2:
    theseclist = [h.a1_111, h.a10_11]
if numbranches==8:
    theseclist = [h.a1_111, h.a8_11, h.a9_122, h.a3_11, h.a10_11, h.a7_1221, h.a5_1, h.a4_121]

sl2 = h.SectionList()
for sec in theseclist:
    sl2.append(sec = sec)
poppedsecs = sl2.unique()
h.refreshnseg(h.makeactivelist(sl2))
print "nseg", h.nsegcnt()
h.cvode.cache_efficient(1)

def getfi(amp,simiter):
    h.tstop = 500 # to match synaptic input runs

    # 'background' injection
    npy.random.seed(int(simiter+h.luckyoffset))
    icr = h.IClamp(h.soma(.5))
    icr.dur = h.tstop
    if numbranches==2:
        icmean = .75
    if numbranches==8:
        icmean = 0
    icstd = 1
    icrand = h.Vector(icmean+icstd*npy.random.randn(h.tstop/h.dt+1))
    icrand.play(icr._ref_amp,h.dt)

    # do current injection run
    ic = h.IClamp(h.soma(.5))
    ap = h.APCount(h.soma(.5))
    ic.delay = 0
    ic.dur = h.tstop 
    ic.amp = amp
    h.run()
    return ap.n*1000./h.tstop

# current amps
imin = 0
if numbranches==2:
    imax = 1.1
if numbranches==8:
    imax = 4
istep = .05
Is = npy.arange(imin,imax+istep,istep)
f = npy.zeros((Is.size))
for ind,I in enumerate(Is):
    f[ind] = getfi(I,simiter)

cPickle.dump(f,open('data/fI-%dbranch-run%d.pkl' % (numbranches,simiter),'w'))
cPickle.dump(Is,open('data/fI-I-%dbranch.pkl' % (numbranches),'w'))
