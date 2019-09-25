#!/usr/bin/env python
from gtools import *
import resource
from datetime import datetime
#_, hard = resource.getrlimit(resource.RLIMIT_DATA)
#resource.setrlimit(resource.RLIMIT_DATA, (2**32,hard))
#print('Max memory usage restricted to 4 GB')

startTime = datetime.now()

names = ['./catalogs/1130788624_source_array.sav']#, './catalogs/1130784064_source_array.sav', './catalogs/1130781304_source_array.sav']
nside = 4096
variances = [0.05, 0.1, 0.15]#, 0.1, 0.03, 0.01]
caps = [1.0, 1.5, 2.0] #[2.0, 1.5, 1.0]
radii = [5.0, 10.0, 15.0, 20.0]
sources = [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1]]
for fname in names:
    for var in variances:
        for cap in caps:
            for radius in radii:
                for src in sources:
                    makeGaussPlot(fname, nside, var, cap, saveData=True, loadData = False, points=src[0], centroids=src[1], components=src[2], cmaps=['Reds','Blues','Greens'], alpha=[1.0,1.0,1.0],radius=radius)
print('tooltest.py executed in '+str(datetime.now()-startTime)+' seconds.')
