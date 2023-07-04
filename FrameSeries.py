#!/usr/bin/env python

'''Runs a series of the Pixel tracking pipelines automatically'''

print("\n*** Pixel Tracking Pipeline Queueing Utility ***\n")

import numpy as np
import sys
import os

dem_file = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem" 

frames = np.loadtxt("frames.txt", dtype=str)

if len(frames) < 2:
    print("Insufficient frames provided")
    sys.exit()
else:
    print(f"{len(frames)} frames provided to compute, {len(frames) - 1} pair(s)y")
    for i in range(len(frames) - 1):
        print(frames[i]+" -> "+frames[i + 1])

response = input('\nRun the pipeline for this series of frames? [y, n]:')
print(response)

for i in range(len(frames) - 1):
        command = f"./PixelTrackingERS.py -f1 {frames[i]} -f2 {frames[i + 1]} -d {dem_file} -dest ~/PixelTrackingRuns/RUN20230704 --ignore-ampcor"

        print("\nRun #"+str(i)+":\n")
        os.system(command)