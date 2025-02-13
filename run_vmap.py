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
parser.add_argument(
    "--convert_disp",
    default=True,
    action=argparse.BooleanOptionalAction,
    dest="convert_disp",
    required=False,
    help="Determines whether to convert the displacements",
)
parser.add_argument(
    "--reproj",
    default=False,
    action=argparse.BooleanOptionalAction,
    dest="reproj",
    required=False,
    help="Determines whether to reproject the output rasters",
)
parser.add_argument(
    "--generate-csv",
    default=True,
    action=argparse.BooleanOptionalAction,
    dest="generateCSV",
    required=False,
    help="Determines whether to convert output rasters to a csv file (point cloud) based on the DEM coordinates.",
)


inps = parser.parse_args()


def generateGeotiff(
    array,
    out_filename,
    folder,
    geometry,
    band_description=None,
    downsample=None,
    **kwargs,
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

    EPSG = 4326  # 8086 #9039  # 3057

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)

    # print(srs.ExportToWkt())

    out_ds.SetProjection(srs.ExportToWkt())

    band = out_ds.GetRasterBand(1)
    band.WriteArray(array)
    if band_description:
        band.SetDescription(band_description)
    band.FlushCache()
    band.ComputeStatistics(False)

    out_ds = None

    os.system(
        f"gdal_translate -r bilinear \
                            -gcp 0 0 {geometry['UpperLeft'][1]} {geometry['UpperLeft'][0]} \
                            -gcp {width} 0 {geometry['UpperRight'][1]} {geometry['UpperRight'][0]} \
                            -gcp {width} {height} {geometry['LowerRight'][1]} {geometry['LowerRight'][0]} \
                            -gcp 0 {height} {geometry['LowerLeft'][1]} {geometry['LowerLeft'][0]} \
                            {os.path.join(folder, out_filename + '_untranslated.tif')} \
                            {os.path.join(folder, out_filename + '_unwarped.tif')}"
    )

    if downsample == True:
        out_name = os.path.join(folder, out_filename + "_uncompressed.tif")
    else:
        out_name = os.path.join(folder, out_filename + ".tif")

    os.system(
        f"gdalwarp \
                   -r bilinear -t_srs EPSG:4326 -et 0 -dstnodata nan \
                  {os.path.join(folder, out_filename + '_unwarped.tif')} \
                  {out_name}"
    )

    os.system(
        f"gdal_translate -co TILED=yes -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 {out_name} {os.path.join(folder, out_filename + '_tiled.tif')}"
    )

    if downsample:
        print("\nExtra downsampling step")

        os.system(
            f"gdal_translate -r bilinear \
                -outsize 3840 0 \
                -co COMPRESS=LERC_ZSTD \
                {os.path.join(folder, out_filename + '_uncompressed.tif')} \
                {os.path.join(folder, out_filename + '.tif')}"
        )

    os.system(f"rm -rf {os.path.join(folder, out_filename + '_untranslated.tif')}")
    os.system(f"rm -rf {os.path.join(folder, out_filename + '_unwarped.tif')}")
    if downsample:
        os.system(f"rm -rf {os.path.join(folder, out_filename + '_uncompressed.tif')}")


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

    print("\n - Reference input file:\n")

    in_ds = gdal.Open(fn1, gdal.GA_ReadOnly)
    in_array = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
    in_ds = None

    out_filename1 = "reference"

    generateGeotiff(in_array, out_filename1, folder1, Geometry)

    ### Secondary

    print("\n - Secondary input file:\n")

    in_ds = gdal.Open(fn2, gdal.GA_ReadOnly)
    in_array = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
    in_ds = None

    out_filename2 = "secondary"

    generateGeotiff(in_array, out_filename2, folder2, Geometry)

    ### Geometry rasters

    print("\n - Generating geometry files in same projection:\n")

    geometry_folder = "geometry_tif"

    os.system(f"rm -rf {geometry_folder}")
    os.mkdir(geometry_folder)

    ### LON

    print(" - Longitude file:\n")

    lon_fn = "geometry/lon.rdr.full"

    in_ds = gdal.Open(lon_fn, gdal.GA_ReadOnly)
    in_array = np.fliplr(in_ds.GetRasterBand(1).ReadAsArray())
    in_ds = None

    out_filename_lon = "lon"

    generateGeotiff(in_array, out_filename_lon, geometry_folder, Geometry)

    ### LAT

    print("\n - Latitude file:\n")

    lat_fn = "geometry/lat.rdr.full"

    in_ds = gdal.Open(lat_fn, gdal.GA_ReadOnly)
    in_array = np.fliplr(in_ds.GetRasterBand(1).ReadAsArray())
    in_ds = None

    out_filename_lat = "lat"

    generateGeotiff(in_array, out_filename_lat, geometry_folder, Geometry)

    ### LAT

    print("\n - Elevation file:\n")

    lat_fn = "geometry/z.rdr.full"

    in_ds = gdal.Open(lat_fn, gdal.GA_ReadOnly)
    in_array = np.fliplr(in_ds.GetRasterBand(1).ReadAsArray())
    in_ds = None

    out_filename_z = "z"

    generateGeotiff(in_array, out_filename_z, geometry_folder, Geometry)

    ### LOS

    print("\n - LOS files:\n")

    lat_fn = "geometry/los.rdr.full"

    in_ds = gdal.Open(lat_fn, gdal.GA_ReadOnly)

    in_array = np.fliplr(in_ds.GetRasterBand(1).ReadAsArray())

    out_filename_inc = "inc"

    generateGeotiff(
        in_array,
        out_filename_inc,
        geometry_folder,
        Geometry,
        band_description="Incidence Angle (+vertical)",
    )

    in_array = np.fliplr(in_ds.GetRasterBand(2).ReadAsArray())

    out_filename_az = "az"

    generateGeotiff(
        in_array,
        out_filename_az,
        geometry_folder,
        Geometry,
        band_description="Azimuth Angle (degrees from North, anti-clockwise)",
    )

    in_ds = None


### location: /home/yadevries/anaconda3/lib/python3.10/site-packages/vmap added to path, make sure to they are executable
## to try: gdal_translate -co TILED=yes -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 input.tif output.tif

if inps.vmap:
    cwd = os.getcwd()

    command = f"vmap.py {os.path.join(cwd, folder1, 'reference.tif')} {os.path.join(cwd, folder2, 'secondary.tif')} -kernel 35 -erode 512 -refinement 2 -dt none"

    print(command)

    os.system(command)
