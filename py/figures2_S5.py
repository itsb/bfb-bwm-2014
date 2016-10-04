# coding: utf-8
import numpy as npy
from time import time
import scipy.interpolate
from mpl_toolkits.mplot3d import axes3d
import cPickle

# set default font sizes, important for figure export
pyl.rc('xtick',labelsize=10)
pyl.rc('ytick',labelsize=10)
pyl.rc('axes',labelsize=10,titlesize=14)
pyl.rc('legend',fontsize=8)
pyl.rc('path',simplify=True)     # http://matplotlib.sourceforge.net/examples/pylab_examples/simplification_clipping_test.html
pyl.rc('figure.subplot',top=.95,bottom=.15,left=.05,right=.95,wspace=.0,hspace=.25)

def computefit(spkr,Itof,ftoI):
    spkrorig = spkr.copy()
    # 2 input (group) regression
    ii, jj = npy.where(spkr>0)

    A = npy.zeros((ii.size,npy.sum(spkr.shape)+1))
    for ind,(i,j) in enumerate(zip(ii,jj+npy.sum(spkr.shape[:1]))):
        A[ind,(i,j)] = 1

    spkr = spkrorig

    # I0 coefficient
    A[:,-1] = 1
    # remove 0 syn coefficients so they are properly absorbed into the I0 coefficient
    A = A[:,list(set(range(A.shape[1]))-set([0,spkr.shape[0]]))]

    b = ftoI(spkr[ii,jj])
    x = npy.linalg.lstsq(A,b)[0]
    best = npy.dot(A,x)
    lmin = npy.min((best,b))
    lmax = npy.max((best,b))

    Ii = npy.array([0]+list(x[:spkr.shape[0]-1]))
    Ij = npy.array([0]+list(x[spkr.shape[0]-1:-1]))
    imin = x[-1]

    iest = npy.zeros(spkr.shape)
    for (i,j) in zip(ii,jj):
        iest[i,j] = Ii[i]+Ij[j]+imin
    spkrest = Itof(iest)

    lmin = npy.min((spkrest,spkr))
    lmax = npy.max((spkrest,spkr))

    mse = ((spkr-spkrest)**2).mean()
    rmse = mse**.5
    nrmse = 100*rmse/spkr.ptp()
    spkcorr = npy.corrcoef(spkr.ravel(),spkrest.ravel())[0,1]
    return (Ii,Ij,imin,spkrest,lmin,lmax,mse,rmse,nrmse,spkcorr)

def f2S5plot(loc1,loc2,numbranches):
    ss = 2 # stepsize
    I = cPickle.load(open('data/fI-I-%dbranch.pkl' % (numbranches)))
    f = cPickle.load(open('data/fI-%dbranch-all.pkl' % (numbranches))).mean(0)
    fn = 'data/paper-%dbranch-%dum_%dum_spkavg.h5' % (numbranches,loc1,loc2)
    if numbranches==2:
        maxsyn2=21
    if numbranches==8:
        maxsyn2=41
    spkr = h5get(fn,'spikes')[:41:ss,:maxsyn2:ss]*5001./5000.

    Itof = scipy.interpolate.interp1d(I,f,kind='linear',bounds_error=False,fill_value=0)
    ftoI = scipy.interpolate.interp1d(f,I,kind='linear',bounds_error=False,fill_value=0)

    Ii,Ij,imin,spkrest,lmin,lmax,mse,rmse,nrmse,spkcorr = computefit(spkr,Itof,ftoI)

    # regression finished, plotting below
    syn1= npy.arange(0,Ij.size*ss,ss)
    syn2= npy.arange(0,Ii.size*ss,ss)
    syns2, syns1 = npy.meshgrid(syn1,syn2)

    h = newfig(figsize=(10,10))    

    ax = h.add_subplot(335, projection='3d')
    ax.plot_surface(syns2,syns1,spkr,cstride=1,rstride=1,cmap=pyl.cm.jet,lw=0)
    ax.view_init(azim=-150,elev=15)
    ax.set_xlabel('Branch A input\n(# of synapses)')
    ax.set_ylabel('Branch B input\n(# of synapses)')
    ax.set_title('Actual Response')
    ax.set_zlim3d((-2,lmax))
    setaspectsquare(ax)

    ax = h.add_subplot(334, projection='3d')
    ax.plot_surface(syns2,syns1,spkrest,cstride=1,rstride=1,cmap=pyl.cm.jet,lw=0)
    ax.view_init(azim=-150,elev=15)
    ax.set_xlabel('Branch A input\n(# of synapses)')
    ax.set_ylabel('Branch B input\n(# of synapses)')
    ax.set_zlabel('Response (Hz)')
    ax.set_title('Predicted Response')
    ax.set_zlim3d((-2,lmax))
    setaspectsquare(ax)

    ax = h.add_subplot(336, projection='3d')
    ax.plot_surface(syns2,syns1,spkr-spkrest,cstride=1,rstride=1,cmap=pyl.cm.jet,lw=0)
    ax.view_init(azim=-150,elev=15)
    ax.set_xlabel('Branch A input\n(# of synapses)')
    ax.set_ylabel('Branch B input\n(# of synapses)')
    ax.set_title('Prediction Errors')
    ax.set_zlim3d((-2,lmax))
    setaspectsquare(ax)

    # current estimates
    # fI curve
    # scatter plots
    # prediction and estimate curves
    # prediction errors

    isubset= I[f<lmax]
    # current model components, current estimates and fI curve
    ax = h.add_subplot(333)
    ax.plot(syn1,Ij,lw=2,c='r',label='Branch A dist stim')
    ax.plot(syn2,Ii,lw=2,c='b',label='Branch B prox stim')
    ax.legend(loc='best')
    ax.set_xlabel('Input (# of synapses)')
    ax.set_ylabel('Estimated current\nreaching soma (nA)')
    ax.set_title('Dendritic i/o curves')
    setaspectsquare(ax)

    ax = h.add_subplot(332)
    ax.plot(isubset,Itof(isubset),lw=2,c='k',ls='-',zorder=2,label='measured fI')
    ax.set_ylim([0,ax.get_ylim()[1]])    
    ax.set_xlabel('Current input (nA)')
    ax.set_ylabel('Response (Hz)')
    ax.set_title('Somatic f-I curve')
    setaspectsquare(ax)

    # prediction vs actual scatter plot
    ax = h.add_subplot(339)

    ms = 10 # marker size for scatter plot
    def sparsify(xdata, ydata,ax = ax,ms=ms,pt=1/72.,histbinfactor=5):
        """ reduces number of data points in scatter plot by binning
            useful for vector graphics export, otherwise optional 
        """
        insize = ax.figure.get_size_inches()*ax.get_window_extent()._get_size()/ax.figure.get_window_extent()._get_size()
        dataspan = npy.max([xdata.ptp(),ydata.ptp()])
        dataperpt = (dataspan/insize).min()*pt
        histbinsize = npy.floor(npy.max([1,dataperpt*ms]))/float(histbinfactor)
        maxedge = npy.max([xdata.max(),ydata.max()])
        minedge = 0
        bins = npy.linspace(minedge,maxedge,(maxedge-minedge)/histbinsize+1)
        hret = npy.histogram2d(xdata.ravel(),ydata.ravel(),bins)[0]
        binsx,binsy = npy.where(hret>0)
        return (bins[binsx], bins[binsy])

    spkrest_sparse, spkr_sparse = sparsify(spkrest,spkr,ax)

    ax.scatter(spkrest_sparse,spkr_sparse,s=ms,c='c',edgecolor='None',zorder=3)
    err = (spkr-spkrest).ravel()
    ax.text(lmax*.05,lmax*.9,'Mean absolute error=%.2f Hz' % (abs(err).mean()),color='k',size=8)

    ax.plot([lmin,lmax],[lmin,lmax],c='.5',dashes=(10,5),lw=3,zorder=4)
    ax.axis((lmin,lmax,lmin,lmax))
    ax.set_xlabel('Predicted response (Hz)')
    ax.set_ylabel('Actual response (Hz)')
    setaspectsquare(ax)

    # companion pdf to go with scatter
    ax = h.add_subplot(331)
    err = (spkr-spkrest).ravel()
    print "err min: %g, max: %g, mean: %g, std: %g" % (err.min(), err.max(), abs(err).mean(), abs(err).std())
    bins = npy.arange(npy.floor(err.min())-.5,npy.ceil(err.max())+.5)
    bincenters = (bins+npy.diff(bins).mean()/2)[:-1]
    normed = True
    p = npy.histogram(err,bins=bins,normed=normed)
    ax.plot(bincenters,p[0],color='c',zorder=1,lw=2)
    ax.set_ylabel('Probability')
    ax.set_xlabel('Prediction Error (Hz)')
    ax.set_xlim([-10,10])
    setaspectsquare(ax)

    # prediction and actual curves
    ax1 = h.add_subplot(337)
    linesp = ax1.plot(syn1,spkrest.T,c='r',ls='--',lw=1.5)
    linesa = ax1.plot(syn1,spkr.T,c='r',lw=2)
    ax1.set_xlabel('Branch A input\n(# of synapses)')
    ax1.set_ylabel('Response (Hz)')
    ax1.set_ylim([0,ax1.get_ylim()[1]])    
    setaspectsquare(ax1)

    ax2 = h.add_subplot(338)
    linesp = ax2.plot(syn2,spkrest,c='b',ls='--',lw=1.5)
    linesa = ax2.plot(syn2,spkr,c='b',lw=2)
    ax2.set_xlabel('Branch B input\n(# of synapses)')
    ax2.set_ylabel('Response (Hz)')
    ax2.legend((linesp[0], linesa[0]), ('Predicted response', 'Actual response'), loc='upper left')
    ax2.set_ylim([0,ax2.get_ylim()[1]])    
    setaspectsquare(ax2)

    if numbranches==2:
        h.savefig('figs/Figure2.png')
    if numbranches==8:
        h.savefig('figs/FigureS5.png')
    return (spkcorr,spkr,spkrest,rmse,nrmse)
