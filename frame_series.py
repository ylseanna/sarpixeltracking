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
    required=False,
    help="Folder to save results to.",
)
inps = parser.parse_args()

dem_file = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem"

comment, frames1, frames2 = np.loadtxt("frames.txt", dtype=str, delimiter=",", unpack=True)

if len(frames1) < 2:
    print("Insufficient frames provided")
    sys.exit()
else:
    print(f"{len(frames1)} pairs provided to compute")
    for i in range(len(frames1)):
        print(frames1[i] + " -> " + frames2[i])

response = input("\nRun the pipeline for this series of frames? [y, n]:")

if response == "y":
    print("Starting series")
    for i in range(len(frames1)): 
        
        command = f"./pixel_tracking_pipeline.py  -f1 {frames1[i]} -f2 {frames2[i]} -d {dem_file} -dest ../vmappixeltrackingresults -pre {comment[i]}_"
        print(command)
        
        print("\nRun #" + str(i) + ":\n")
        os.system(command)
else:
    print("Aborted")
