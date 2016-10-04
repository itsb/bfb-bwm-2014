from time import time
import numpy as npy
from socket import gethostname

# set default font sizes, important for figure export
pyl.rc('xtick',labelsize=10)
pyl.rc('ytick',labelsize=10)
pyl.rc('axes',labelsize=10,titlesize=14)
pyl.rc('legend',fontsize=8)
pyl.rc('path',simplify=True)     # http://matplotlib.sourceforge.net/examples/pylab_examples/simplification_clipping_test.html
pyl.rc('figure.subplot',top=.95,bottom=.15,left=.05,right=.95,wspace=.0,hspace=.25)

from neuron import h
h.load_file('basal_project.hoc')

def postrunrecgather(vd):
    ags = []
    for a in h.ampaglist:
        ags.append(npy.array(a))
    ags = npy.array(ags)
    vd.update({'ags':ags})

    ngs = []
    for n in h.nmdaglist:
        ngs.append(npy.array(n))
    ngs = npy.array(ngs)
    vd.update({'ngs':ngs})

    ais = []
    for a in h.ampailist:
        ais.append(npy.array(a))
    ais = npy.array(ais)
    vd.update({'ais':ais})

    nis = []
    for n in h.nmdailist:
        nis.append(npy.array(n))
    nis = npy.array(nis)
    vd.update({'nis':nis})

    nbs = []
    for n in h.nmdablist:
        nbs.append(npy.array(n))
    nbs = npy.array(nbs)
    vd.update({'nbs':nbs})

    rasind = []
    for s in h.ex_stim_vecs:
        rasind.append(npy.array(s))
    vd.update({'rasind':rasind})

def pyrun(ratio=0, loc=0, nsyn=0, sec=h.a10_11, icamps=[],dendrec=True):

    print '\n%s starting run' % (gethostname())
    iotim = 0
    tic = time()

    h.tstop = 500

    h.tsamp = h.tstop/h.dt+1
    h.synrec = 1
    tsamp = int(h.tsamp)

    r = h.Random(h.luckyoffset)
    r.negexp(h.meanisi)

    # stimulated branch
    syn = h.synsuper(.5, r, sec=sec)
    syndict = dict(sloc=loc, ratio=ratio, e1flag=0)
    for name, value in syndict.iteritems(): 
        setattr(syn.syn, name, value)

    # initialize nseg with 'active' branches
    seclist = [h.a1_111, h.a10_11]
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
        d = h.distance(0, sec=sec)
        locx = (loc - d)/sec.L
        v = h.Vector(tsamp)
        v.record(sec(locx)._ref_v)
        trash = v.label(sec.name())
        vd.update({'d':v})

    h.poisson = 1
    # 'background' current injection
    ic = h.IClamp(0.5,sec=h.soma)
    ic.dur = h.tstop

    p_tstart = 100 # plot tstart, crop initial rising phase of current injection
    ind_tstart = p_tstart/h.dt
    t = npy.arange(0,h.tstop-p_tstart+h.dt,h.dt)
    fh = newfig(figsize=(4,8))
    for runcnt, icamp in enumerate(icamps):
        ic.amp = icamp
        syn.syn.nsyn = nsyn
        seed1 = float(686)
        r1 = h.Random(seed1)
        r1.negexp(h.meanisi)
        syn.setrand(r1)

        # run simulation
        h.run()

        postrunrecgather(vd)

        # plot voltage and synaptic current
        ax = fh.add_subplot(6,1,runcnt+3)
        ax.plot(t,npy.array(vd['d'])[ind_tstart:],c='k',lw=1)
        ax = fh.add_subplot(6,1,runcnt+5)
        ax.plot(t,-(vd['nis'].sum(0)+vd['ais'].sum(0))[ind_tstart:],c='r',lw=1)

    ax = fh.add_subplot(6,1,1)
    for ind,ras in enumerate(vd['rasind']):
        ticks = ras[ras>=p_tstart]
        ax.scatter(ticks-p_tstart,npy.ones(ticks.shape)*ind,marker=[[[0,0],[0,1]],0])
    ax.axis((0,h.tstop-p_tstart,0,ind+1))

    zero = npy.zeros(t.size)
    ax = fh.add_subplot(6,1,2)
    ax.fill_between(t,zero,(vd['ngs']/vd['nbs']).sum(0)[ind_tstart:],lw=2,color='k')
    ax.fill_between(t,zero,vd['ags'].sum(0)[ind_tstart:],lw=2,color='r')
    ax.set_xlim([0,h.tstop-p_tstart])

    fh.savefig('figs/Figure4.png')

    print '%s running %d runs took %d seconds' % (
           gethostname(), runcnt+1, time() - tic)
