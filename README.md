# SAR Pixel Tracking Pipeline

Mirror of an original SAR pixel tracking pipeline for a Rannís funded project (located [here](https://github.com/ylseanna/RANNIS2023PixelTracking)).

A number of scripts to run AutoRIFT with ISCE for PixelTracking applications on the University of Iceland's Dyngja computer

The database scripts can be adjusted to create a small sqlite database of files to be processed, this is then used by the main Pixel Tracking Script

Currently it is set up to just ingest data from a particular folder

## To-do

- Make file inputs a command line argument
- Give command line argument for doing DenseOffsets
- Give Dem input otherwise use srtm 
