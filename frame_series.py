#!/usr/bin/env python

"""Runs a series of the Pixel tracking pipelines automatically"""

print("\n*** Pixel Tracking Pipeline Queueing Utility ***\n")

import numpy as np
import sys
import os

import argparse

parser = argparse.ArgumentParser(
    prog="ERS Pixel Tracking Pipeline Queueing Utility",
    description="Completes a full run of the ERS pixel tracking pipelines for a series of frames (stripmap) images.",
)
parser.add_argument(
    "-f",
    "--folder",
    dest="folder",
    type=str,
    required=True,
    help="Folder to save results to.",
)
inps = parser.parse_args()

dem_file = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem"

frames = np.loadtxt("frames.txt", dtype=str)

if len(frames) < 2:
    print("Insufficient frames provided")
    sys.exit()
else:
    print(f"{len(frames)} frames provided to compute, {len(frames) - 1} pair(s)")
    for i in range(len(frames) - 1):
        print(frames[i] + " -> " + frames[i + 1])

response = input("\nRun the pipeline for this series of frames? [y, n]:")

if response == "y":
    print("Starting series")
    for i in range(len(frames) - 1):
        command = f"./pixel_tracking.py -f1 {frames[i]} -f2 {frames[i + 1]} -d {dem_file} -dest ~/PixelTrackingRuns/{inps.folder} --ignore-ampcor"

        print("\nRun #" + str(i) + ":\n")
        os.system(command)
else:
    print("Aborted")
