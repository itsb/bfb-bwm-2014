#!/usr/bin/env python

# suprathreshold simulation config
nmpa = 2.6 # NMDA-AMPA ratio
loc = 90
nsyn = 20
icamps = [1.05,1.73]

# launch simulations
execfile('py/common.py')
execfile('hoc/pyrun.py')
pyrun(ratio=nmpa,loc=loc,nsyn=nsyn,icamps=icamps)
