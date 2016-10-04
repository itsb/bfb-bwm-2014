from time import time
import os.path
import tables as pyt
import numpy as npy
from socket import gethostname

from neuron import h
h.load_file('basal_project.hoc')

def pyloop(ratio=0, loc1=128, loc2=128, simiter=0, numbranches=0, dendrec=False):

    print '\n%s starting run' % (gethostname())
    iotim = 0
    tic = time()

    h5fn = 'paper-%dbranch-%dum_%dum_run%d.h5' % (numbranches, l1, l2, simiter)
    h.tstop = 500

    savd = './data/'
    filters = pyt.Filters(complevel=9, complib='zlib', shuffle=True, fletcher32=True)
    h5f = pyt.openFile(os.path.join(savd, h5fn),'w',filters=filters)

    maxsyn1 = int(h.maxsyn1)+1
    if numbranches==2:
        h.maxsyn2 = 20
    if numbranches==8:
        h.maxsyn2 = 40
    maxsyn2 = int(h.maxsyn2)+1
    h.tsamp = h.tstop/h.dt+1
    tsamp = int(h.tsamp)
    shape = (maxsyn1, maxsyn2, tsamp)
    cshape = (1, 1, shape[-1])

    r = h.Random(simiter+h.luckyoffset)
    r.negexp(h.meanisi)

    # stimulated branch
    if numbranches==2:
        seclist = [[h.a1_111], [h.a10_11]]
    if numbranches==8:
        seclist = [[h.a1_111, h.a8_11, h.a9_122, h.a3_11],[h.a10_11, h.a7_1221, h.a5_1, h.a4_121]]

    syns = [[], []]
    for ind, loc in enumerate([loc1,loc2]):
        for sec in seclist[ind]:
            syn = h.synsuper(.5, r, sec=sec)
            syns[ind].append(syn)
            syndict = dict(sloc=loc, ratio=ratio, e1flag=0)
            for name, value in syndict.iteritems(): 
                setattr(syn.syn, name, value)

    # initialize nseg with 'active' branches
    seclist = [sec for sublist in seclist for sec in sublist]
    sl2 = h.SectionList()
    for sec in seclist:
        sl2.append(sec = sec)
    poppedsecs = sl2.unique()
    h.refreshnseg(h.makeactivelist(sl2))
    print 'nseg: %d' % (h.nsegcnt())
    h.cvode.cache_efficient(1)
    h.cvode_active(0)

    # somatic voltage recordings, after nseg initialization
    v = h.Vector(tsamp)
    v.record(h.soma(.5)._ref_v)
    trash = v.label(h.soma.name())
    # voltage recording dictionary
    vd = {'s':v}

    if dendrec:
        # dendritic voltage recording 
        h.distance(0, h.soma_con_pt, sec=h.soma)
        for n, sec in enumerate(seclist):
            d = h.distance(0, sec=sec)
            # location 1 synapses
            locx1 = (loc1 - d)/sec.L
            v = h.Vector(tsamp)
            v.record(sec(locx1)._ref_v)
            trash = v.label(sec.name())
            vd.update({'dp'+str(n):v})
            # location 2 synapses
            locx2 = (loc2 - d)/sec.L
            v = h.Vector(tsamp)
            v.record(sec(locx2)._ref_v)
            trash = v.label(sec.name())
            vd.update({'dd'+str(n):v})

    h.poisson = 1
    # 'background' current injection
    npy.random.seed(int(simiter+h.luckyoffset))
    ic = h.IClamp(0.5,sec=h.soma)
    ic.dur = h.tstop
    if numbranches==2:
        icmean = .75
    if numbranches==8:
        icmean = 0
    icstd = 1
    icvals = icmean+icstd*npy.random.randn(h.tstop/h.dt+1)
    icrand = h.Vector(tsamp)
    for i in xrange(tsamp):
        icrand.x[i] = icvals[i]
    icrand.play(ic._ref_amp, h.dt)

    h5d = {}
    for key in vd.keys():
        h5d.update({key:h5f.createCArray(h5f.root, key, shape=shape,
                                         atom=pyt.Float64Atom(),
                                         title=vd[key].label(),
                                         chunkshape=cshape)
                                         })

    synpairs = [(i, j) for i in xrange(0,maxsyn1,2) for j in xrange(0,maxsyn2,2)]
    iisi = simiter*(maxsyn1 + maxsyn2)  # inter-iteration seed interval

    for runcnt, (nsyn1, nsyn2) in enumerate(synpairs):
        # configure s(t)imulation
        for ssyn in syns[0]: 
            ssyn.syn.nsyn = nsyn1 
            seed1 = float(iisi + nsyn1 + h.luckyoffset)
            r1 = h.Random(seed1)
            r1.negexp(h.meanisi)
            ssyn.setrand(r1)
        for ssyn in syns[1]: 
            ssyn.syn.nsyn = nsyn2 
            seed2 = float(iisi + maxsyn1 + nsyn2 + h.luckyoffset)
            r2 = h.Random(seed2)
            r2.negexp(h.meanisi)
            ssyn.setrand(r2)

        # run simulation
        h.run()

        # write results
        iotic = time()
        for key in h5d.keys():
            h5d[key][nsyn1, nsyn2] = npy.array(vd[key]).reshape(cshape)
        iotim += time() - iotic

    iotic = time()
    h5f.close()
    iotim += time() - iotic

    print '%s running %d runs took %d seconds, of which %d seconds was I/O' % (
           gethostname(), runcnt+1, time() - tic, iotim)
