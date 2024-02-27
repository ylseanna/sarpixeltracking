#!/usr/bin/env python

import os
from osgeo import gdal
import numpy as np

import argparse

parser = argparse.ArgumentParser(
    prog="vmap ISCE compat",
    description="runs vmap with ISCE",
)
parser.add_argument(
    "-f1",
    "--file1",
    dest="file1",
    type=str,
    required=False,
    help="Input filename 1",
)
parser.add_argument(
    "-f2",
    "--file2",
    dest="file2",
    type=str,
    required=False,
    help="Input filename 2",
)
parser.add_argument(
    "--geocode",
    default=True,
    action=argparse.BooleanOptionalAction,
    dest="geocode",
    required=False,
    help="Determines whether to generate geocoded tiffs",
)
parser.add_argument(
    "--vmap",
    default=True,
    action=argparse.BooleanOptionalAction,
    dest="vmap",
    required=False,
    help="Determines whether to run vmap",
)

inps = parser.parse_args()


def generateGeotiffvmap(
    array, out_filename, folder, geometry, **kwargs
):  # no extension on filename
    import os
    from osgeo import gdal, osr
    import numpy as np

    driver = gdal.GetDriverByName("GTiff")

    out_ds = driver.Create(
        os.path.join(folder, out_filename + "_untranslated.tif"),
        array.shape[1],
        array.shape[0],
        1,
        gdal.GDT_Float32,
    )

    height, width = np.shape(array)

    EPSG = 3857 #8086 #9039  # 3057
    
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)

    print(srs.ExportToWkt())

    out_ds.SetProjection(srs.ExportToWkt())

    band = out_ds.GetRasterBand(1)
    band.WriteArray(array)
    band.FlushCache()
    band.ComputeStatistics(False)

    out_ds = None

    os.system(
        f"gdal_translate -r bilinear \
                            -gcp 0 0 {geometry['UpperLeft'][1]} {geometry['UpperLeft'][0]} \
                            -gcp {width} 0 {geometry['UpperRight'][1]} {geometry['UpperRight'][0]} \
                            -gcp {width} {height} {geometry['LowerRight'][1]} {geometry['LowerRight'][0]} \
                            -gcp 0 {height} {geometry['LowerLeft'][1]} {geometry['LowerLeft'][0]} \
                            {os.path.join(folder, out_filename+'_untranslated.tif')} \
                            {os.path.join(folder, out_filename+'_unwarped.tif')}"
    )

    os.system(
        f"gdalwarp \
                   -r bilinear -t_srs EPSG:{EPSG} -et 0 -dstnodata nan \
                  {os.path.join(folder, out_filename+'_unwarped.tif')} \
                  {os.path.join(folder, out_filename+'_untiled.tif')}"
    )

    os.system(
        f"gdal_translate -co TILED=yes -co BLOCKXSIZE=128 -co BLOCKYSIZE=256 {os.path.join(folder, out_filename+'_untiled.tif')} {os.path.join(folder, out_filename+'.tif')}"
    )


    # os.system(f"rm -rf {os.path.join(folder, out_filename+'_untranslated.tif')}")
    # os.system(f"rm -rf {os.path.join(folder, out_filename+'_unwarped.tif')}")
    # os.system(f"rm -rf {os.path.join(folder, out_filename+'_untiled.tif')}")



def generateGeometry():
    from osgeo import gdal
    import numpy as np
    from geographiclib.geodesic import Geodesic as geodesic

    print("\n - Getting geometric bounds:\n")

    lats = gdal.Open("geometry/lat.rdr.full")
    lons = gdal.Open("geometry/lon.rdr.full")

    lat = np.array(lats.GetRasterBand(1).ReadAsArray())
    lon = np.array(lons.GetRasterBand(1).ReadAsArray())

    lats = None
    lons = None

    height, width = np.shape(lat)

    UpperLeft = (lat[0, -1], lon[0, -1])
    UpperRight = (lat[0, 0], lon[0, 0])
    LowerLeft = (lat[-1, -1], lon[-1, -1])
    LowerRight = (lat[-1, 0], lon[-1, 0])

    print("Upper left corner: ", UpperLeft)
    print("Upper right corner:", UpperRight)
    print("Lower right corner:", LowerRight)
    print("Lower left corner: ", LowerLeft)

    top_geodesic = geodesic.WGS84.Inverse(
        UpperLeft[0], UpperLeft[1], UpperRight[0], UpperRight[1]
    )
    bottom_geodesic = geodesic.WGS84.Inverse(
        LowerLeft[0], LowerLeft[1], LowerRight[0], LowerRight[1]
    )
    left_geodesic = geodesic.WGS84.Inverse(
        UpperLeft[0], UpperLeft[1], LowerLeft[0], LowerLeft[1]
    )
    right_geodesic = geodesic.WGS84.Inverse(
        UpperRight[0], UpperRight[1], LowerRight[0], LowerRight[1]
    )

    print("\nTop geodesic distance:   ", top_geodesic["s12"])
    print("Bottom geodesic distance:", bottom_geodesic["s12"])
    print("Left geodesic distance:  ", left_geodesic["s12"])
    print("Right geodesic distance: ", right_geodesic["s12"])

    Geometry = {
        "UpperLeft": UpperLeft,
        "UpperRight": UpperRight,
        "LowerRight": LowerRight,
        "LowerLeft": LowerLeft,
        "left_geodesic": left_geodesic,
        "right_geodesic": right_geodesic,
        "top_geodesic": top_geodesic,
        "bottom_geodesic": bottom_geodesic,
        "pixel_height": height,
        "pixel_width": width,
    }

    return Geometry


folder1 = "reference_tif"
folder2 = "coreg_secondary_tif"

import json

log = open("log.json")
logdata = json.load(log)

reference = logdata["frame_metadata"][0]["reference"]
secondary = logdata["frame_metadata"][1]["secondary"]

referencetimestring = reference["begintime"][:10].replace("-", "")
secondarytimestring = secondary["begintime"][:10].replace("-", "")

print(referencetimestring, secondarytimestring)

if inps.geocode == True:
    Geometry = generateGeometry()

    os.system(f"rm -rf {folder1}")
    os.mkdir(folder1)
    os.system(f"rm -rf {folder2}")
    os.mkdir(folder2)

    if os.path.exists("reference_slc_crop"):
        fn1 = "reference_slc_crop/reference.slc"
    else:
        fn1 = "reference_slc/reference.slc"

    fn2 = "coregisteredSlc/refined_coreg.slc"

    ### Reference

    in_ds = gdal.Open(fn1, gdal.GA_ReadOnly)
    in_array = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
    in_ds = None

    out_filename1 = "reference_" + referencetimestring

    generateGeotiffvmap(in_array, out_filename1, folder1, Geometry)

    ### Secondary

    in_ds = gdal.Open(fn2, gdal.GA_ReadOnly)
    in_array = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
    in_ds = None

    out_filename2 = "secondary_" + secondarytimestring

    generateGeotiffvmap(in_array, out_filename2, folder2, Geometry)


### location: /home/yadevries/anaconda3/lib/python3.10/site-packages/vmap added to path, make sure to they are executable
## to try: gdal_translate -co TILED=yes -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 input.tif output.tif

if inps.vmap == True:
    cwd = os.getcwd()

    command = f"vmap.py {os.path.join(cwd, folder1, 'reference_'+referencetimestring+'.tif')} {os.path.join(cwd, folder2, 'secondary_'+secondarytimestring+'.tif')} -kernel 35 -erode 1024 -refinement 2"

    print(command)

    os.system(command)
