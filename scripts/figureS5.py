#!/usr/bin/env python
import sys

# suprathreshold simulation config
runnum = int(sys.argv[1])-1
nmpa = 2.6 # NMDA-AMPA ratio
l1, l2 = 90, 150

# launch simulations
execfile('hoc/pyloop.py')
pyloop(ratio=nmpa, loc1=l1, loc2=l2, simiter=runnum, numbranches=8)
