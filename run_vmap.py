#!/usr/bin/env python

import os
from osgeo import gdal, osr
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


def generateGeotiffvmap(
    array, out_filename, folder, geometry, band_description=None, **kwargs
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

    EPSG = 4326 #8086 #9039  # 3057
    
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
                            {os.path.join(folder, out_filename+'_untranslated.tif')} \
                            {os.path.join(folder, out_filename+'_unwarped.tif')}"
    )

    os.system(
        f"gdalwarp \
                   -r bilinear -t_srs EPSG:{EPSG} -et 0 -dstnodata nan \
                  {os.path.join(folder, out_filename+'_unwarped.tif')} \
                  {os.path.join(folder, out_filename+'_untiled.tif')}"
    )
    
    output_EPSG=4326

    os.system(
        f"gdal_translate -co TILED=yes -co BLOCKXSIZE=128 -co BLOCKYSIZE=256 -a_srs EPSG:{output_EPSG} {os.path.join(folder, out_filename+'_untiled.tif')} {os.path.join(folder, out_filename+'.tif')}"
    )


    os.system(f"rm -rf {os.path.join(folder, out_filename+'_untranslated.tif')}")
    os.system(f"rm -rf {os.path.join(folder, out_filename+'_unwarped.tif')}")
    os.system(f"rm -rf {os.path.join(folder, out_filename+'_untiled.tif')}")



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

print("\n - Frame metadata:\n")

print(reference)

print(secondary)



if inps.geocode == True:  
    print("\n - Generating vmap input files (simple initial projection):\n")

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

    generateGeotiffvmap(in_array, out_filename1, folder1, Geometry)

    ### Secondary
    
    print("\n - Secondary input file:\n")


    in_ds = gdal.Open(fn2, gdal.GA_ReadOnly)
    in_array = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
    in_ds = None

    out_filename2 = "secondary"

    generateGeotiffvmap(in_array, out_filename2, folder2, Geometry)
    
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

    generateGeotiffvmap(in_array, out_filename_lon, geometry_folder, Geometry)
    
    ### LAT

    print("\n - Latitude file:\n")

    lat_fn = "geometry/lat.rdr.full"

    in_ds = gdal.Open(lat_fn, gdal.GA_ReadOnly)
    in_array = np.fliplr(in_ds.GetRasterBand(1).ReadAsArray())
    in_ds = None

    out_filename_lat = "lat"

    generateGeotiffvmap(in_array, out_filename_lat, geometry_folder, Geometry)
    
    ### LAT

    print("\n - Elevation file:\n")

    lat_fn = "geometry/z.rdr.full"

    in_ds = gdal.Open(lat_fn, gdal.GA_ReadOnly)
    in_array = np.fliplr(in_ds.GetRasterBand(1).ReadAsArray())
    in_ds = None

    out_filename_z = "z"

    generateGeotiffvmap(in_array, out_filename_z, geometry_folder, Geometry)
    
    ### LOS

    print("\n - LOS files:\n")

    lat_fn = "geometry/los.rdr.full"

    in_ds = gdal.Open(lat_fn, gdal.GA_ReadOnly)
    
    in_array = np.fliplr(in_ds.GetRasterBand(1).ReadAsArray())
    
    out_filename_inc = "inc"
    
    generateGeotiffvmap(in_array, out_filename_inc, geometry_folder, Geometry, band_description="Incidence Angle (+vertical)")
    
    in_array = np.fliplr(in_ds.GetRasterBand(2).ReadAsArray())
    
    out_filename_az = "az"
    
    generateGeotiffvmap(in_array, out_filename_az, geometry_folder, Geometry, band_description="Azimuth Angle (degrees from North, anti-clockwise)")
    
    in_ds = None

    

    
    


### location: /home/yadevries/anaconda3/lib/python3.10/site-packages/vmap added to path, make sure to they are executable
## to try: gdal_translate -co TILED=yes -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 input.tif output.tif

if inps.vmap == True:
    cwd = os.getcwd()

    command = f"vmap.py {os.path.join(cwd, folder1, 'reference.tif')} {os.path.join(cwd, folder2, 'secondary.tif')} -kernel 35 -erode 512 -refinement 2 -dt none"

    print(command)

    os.system(command)
    
    
if inps.convert_disp == True:
    from pyproj import Transformer

    from glob import glob
    from pygeotools.lib import iolib, geolib, warplib, timelib
    from geographiclib.geodesic import Geodesic
    
    print("\n - Converting vmap displacements\n")
    
    vmap_dir = glob("reference*vmap*")[0]
    
    ds = iolib.fn_getds(vmap_dir+"/vmap-F.tif")
    
    dX_pix, dY_pix, goodPixels = ds.ReadAsArray()
    
    ref_ds = iolib.fn_getds("reference_tif/reference.tif")
    
    data = ref_ds.ReadAsArray()
    
    num_pix_y, num_pix_x = data.shape
    
    print(num_pix_x, num_pix_y)
    
    extent = geolib.ds_extent(ref_ds)
    
    print(extent)
    
    UpperLeft = [extent[0], extent[3]]
    UpperRight = [extent[2], extent[3]]
    LowerLeft = [extent[0], extent[1]]
    LowerRight = [extent[2], extent[1]]
    
    
    # top_geodesic = Geodesic.WGS84.Inverse(
    #     UpperLeft[0], UpperLeft[1], UpperRight[0], UpperRight[1]
    # )
    # bottom_geodesic = Geodesic.WGS84.Inverse(
    #     LowerLeft[0], LowerLeft[1], LowerRight[0], LowerRight[1]
    # )
    # left_geodesic = Geodesic.WGS84.Inverse(
    #     UpperLeft[0], UpperLeft[1], LowerLeft[0], LowerLeft[1]
    # )
    # right_geodesic = Geodesic.WGS84.Inverse(
    #     UpperRight[0], UpperRight[1], LowerRight[0], LowerRight[1]
    # )
    
    # print("\nTop geodesic distance:   ", top_geodesic["s12"])
    # print("Bottom geodesic distance:", bottom_geodesic["s12"])
    # print("Left geodesic distance:  ", left_geodesic["s12"])
    # print("Right geodesic distance: ", right_geodesic["s12"])
    
    
    WGS86toISN2016 = Transformer.from_crs("EPSG:4326", "EPSG:8088")
    
    UpperLeftISN = WGS86toISN2016.transform(extent[3], extent[0])
    UpperRightISN = WGS86toISN2016.transform(extent[3], extent[2])
    LowerLeftISN = WGS86toISN2016.transform(extent[1], extent[0])
    LowerRightISN = WGS86toISN2016.transform(extent[1], extent[2])
    
    print(UpperLeftISN, UpperRightISN, LowerLeftISN, LowerRightISN)
    
    print("\nTop geodesic distance:   ",  UpperRightISN[0] - UpperLeftISN[0])
    print("Bottom geodesic distance:", LowerRightISN[0] - LowerLeftISN[0])
    print("Left geodesic distance:  ", UpperLeftISN[1] - LowerLeftISN[1])
    print("Right geodesic distance: ", UpperRightISN[1] - LowerRightISN[1])
    
    TopDist = UpperRightISN[0] - UpperLeftISN[0]
    BottomDist = LowerRightISN[0] - LowerLeftISN[0]
    LeftDist = UpperLeftISN[1] - LowerLeftISN[1]
    RightDist =  UpperRightISN[1] - LowerRightISN[1]
    
        
        
    resX = ((TopDist+BottomDist)/2) / num_pix_x
    
    resY = ((LeftDist + RightDist)/2) / num_pix_y
    
        
    # resX = ((top_geodesic["s12"]+bottom_geodesic["s12"])/2) / num_pix_x
    
    # resY = ((left_geodesic["s12"]+right_geodesic["s12"])/2) / num_pix_y
    
    print(resX)
    print(resY)
    
    dX_real = dX_pix * resX
    dY_real = dY_pix * resY
    
    copy_ds = iolib.fn_getds("reference_tif/reference.tif")
    
    iolib.writeGTiff(dX_real,
                     vmap_dir+"/vmap-dx.tif",
                     src_ds=copy_ds,
                     bnum=1)

    iolib.writeGTiff(dY_real,
                     vmap_dir+"/vmap-dy.tif",
                     src_ds=copy_ds,
                     bnum=1)


if inps.reproj == True:
    print("\n - Reprojection to DEM for plotting (using GCPS):\n")
    
    print(" - Defining GCPs:\n")

    ds = gdal.Open("/home/yad2/sarpixeltracking/geometry_tif/lon.tif")
    lon_array = ds.GetRasterBand(1).ReadAsArray()
    maxlon, minlon = lon_array.max(), lon_array.min()
    ds=None
    
    ds = gdal.Open("/home/yad2/sarpixeltracking/geometry_tif/lat.tif")
    lat_array = ds.GetRasterBand(1).ReadAsArray()
    maxlat, minlat = lat_array.max(), lat_array.min()
    ds=None
    
    GCPs = []
    
    num_lines, num_pixels = np.shape(lon_array)
    
    grid_spacing = 50
    
    for i in range(0,num_lines, grid_spacing):
        for j in range(0, num_pixels, grid_spacing):
            lon = float(lon_array[i,j])
            lat = float(lat_array[i,j])
                        
            if not np.isnan(lon) and not np.isnan(lat):
                print("GCPs:", lon, lat, j, i)
                                
                gcp = gdal.GCP(lon, lat, 0, j, i)
                
                GCPs.append(gcp)

    print("Number of GCPs:", len(GCPs))
    
    print("\nSpatial reference for GCPs:")
    
    EPSG = 4326 #8086 #9039  # 3057
    
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    
    print(srs)
        
    GCP_wkt = srs.ExportToWkt()
    
    srs=None
    
    print("\n - Reprojection of reference:\n")
    
    root = "reference"
    
    reproj_folder = f"{root}_reproj"
    
    os.system(f"rm -rf {reproj_folder}")
    os.mkdir(f"{reproj_folder}")
    
    ds = gdal.Open(f"{root}_tif/{root}_{referencetimestring}.tif")
    ds.SetGCPs(GCPs, GCP_wkt)
    
    dstDS = gdal.Warp(f"{reproj_folder}/{root}_reproj.tif", ds, format='GTiff', tps=True, resampleAlg=gdal.GRA_Bilinear)
    
    ds=None
    dstDS=None
    
    
if inps.generateCSV == True:
    from glob import glob
    import pandas as pd
    from pygeotools.lib import iolib, geolib
    from tqdm import tqdm
    
    print("\n - Generating CSV file based on vmap outputs:\n")
    
    print(" - Importing outputs...\n")

    ds = iolib.fn_getds("geometry_tif/lon.tif")
    lon_array = ds.ReadAsArray()
    ds=None
    
    ds = iolib.fn_getds("geometry_tif/lat.tif")
    lat_array = ds.ReadAsArray()
    ds=None
    
    ds = iolib.fn_getds("geometry_tif/z.tif")
    z_array = ds.ReadAsArray()
    ds=None
    
    ds = iolib.fn_getds("geometry_tif/az.tif")
    az_array = ds.ReadAsArray()
    ds=None
    
    ds = iolib.fn_getds("geometry_tif/inc.tif")
    inc_array = ds.ReadAsArray()
    ds=None
    
    #### vmap outputs
    
    vmap_dir = glob("reference*vmap*")[0]
    
    ds = iolib.fn_getds(vmap_dir+"/vmap-F.tif")
    
    dX_pix_array, dY_pix_array, good_pixels = ds.ReadAsArray()
    
    ds=None
    
    ds = iolib.fn_getds(vmap_dir+"/vmap-dx.tif")
    
    dX_m_array = ds.ReadAsArray()
    
    ds=None
    
    ds = iolib.fn_getds(vmap_dir+"/vmap-dy.tif")
    
    dY_m_array = ds.ReadAsArray()
    
    ds=None
    
    ### Eventually not hard-code
    
    extent = [-19.28, -19.41, 63.637, 63.67]
    
    ### generate table
    
    array_shape = good_pixels.shape
    
    temp = []
    
    print(" - Importing looping over arrays:\n")
    
    # for i in range(array_shape[0]):
    for i in tqdm(range(array_shape[0]), desc="Row"):
        for j in range(array_shape[1]):
            if good_pixels[i,j] == 1:
                # print(lat_array[i,j], extent[3])
                # print(lat_array[i,j], extent[2])
                # print(lon_array[i,j], extent[1])
                # print(lon_array[i,j], extent[0])
                
                if lat_array[i,j] < extent[3] and lat_array[i,j] > extent[2] and lon_array[i,j] > extent[1] and lon_array[i,j] < extent[0]:
                    temp.append({
                        "x_disp_p": dX_pix_array[i,j],
                        "y_disp_p": dY_pix_array[i,j],
                        "x_disp_m": dX_m_array[i,j],
                        "y_disp_m": dY_m_array[i,j],
                        "lat": lat_array[i,j],
                        "lon": lon_array[i,j],
                        "z": z_array[i,j],
                        "az": az_array[i,j],
                        "inc": inc_array[i,j],
                    })
    
    df = pd.DataFrame.from_dict(temp)
    
    print(df)
    
    df.to_csv("vmap-offsets.csv", index=False)
    
    # GCPs = []
    
    # num_lines, num_pixels = np.shape(lon_array)
    
    # grid_spacing = 50
    
    # for i in range(0,num_lines, grid_spacing):
    #     for j in range(0, num_pixels, grid_spacing):
    #         lon = float(lon_array[i,j])
    #         lat = float(lat_array[i,j])
                        
    #         if not np.isnan(lon) and not np.isnan(lat):
    #             print("GCPs:", lon, lat, j, i)
                                
    #             gcp = gdal.GCP(lon, lat, 0, j, i)
                
    #             GCPs.append(gcp)

    # print("Number of GCPs:", len(GCPs))
    
    # print("\nSpatial reference for GCPs:")
    
    # EPSG = 4326 #8086 #9039  # 3057
    
    # srs = osr.SpatialReference()
    # srs.ImportFromEPSG(EPSG)
    
    # print(srs)
        
    # GCP_wkt = srs.ExportToWkt()
    
    # srs=None
    
    # print("\n - Reprojection of reference:\n")
    
    # root = "reference"
    
    # reproj_folder = f"{root}_reproj"
    
    # os.system(f"rm -rf {reproj_folder}")
    # os.mkdir(f"{reproj_folder}")
    
    # ds = gdal.Open(f"{root}_tif/{root}_{referencetimestring}.tif")
    # ds.SetGCPs(GCPs, GCP_wkt)
    
    # dstDS = gdal.Warp(f"{reproj_folder}/{root}_reproj.tif", ds, format='GTiff', tps=True, resampleAlg=gdal.GRA_Bilinear)
    
    # ds=None
    # dstDS=None
    
