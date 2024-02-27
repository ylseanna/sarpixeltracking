#!/usr/bin/env python

import os
from osgeo import gdal, osr
import netCDF4 as nc
import numpy as np

import xml.etree.ElementTree

def generate_simple_GeoTIFF(in_file, out_file, flip_ud=False):
    in_ds = gdal.Open(in_file, gdal.GA_ReadOnly)
    band  = in_ds.GetRasterBand(1)
    array = band.ReadAsArray()
    in_ds = None
    
    cols, rows = array.shape
    
    if flip_ud == True:
        array = np.flipud(array)

    dtype = band.DataType
    
    driver  = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(out_file, rows, cols, 1, dtype)
    outband = outdata.GetRasterBand(1)
    outband.WriteArray(array)
    outdata = None


def generate_geocoded_tif(array, file_name_root):
    os.system(f"rm -rf {file_name_root}_amp")
    os.mkdir(file_name_root + "_amp")

    cols, rows = array.shape

    print(rows, cols)

    print("\n - Generating amplitude file")
    amp_filename = os.path.join(
        os.getcwd(), file_name_root + "_amp", file_name_root + "_amp.tif"
    )
    print(amp_filename)

    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(amp_filename, rows, cols, 3, gdal.GDT_Float32)
    outband = outdata.GetRasterBand(1)
    outband.WriteArray(array)
    
    ds = gdal.Open("/home/yad2/sarpixeltracking/geometry_tif/lon.tif")
    lon_array = ds.GetRasterBand(1).ReadAsArray()
    ds=None
    
        
    lonband = outdata.GetRasterBand(2)
    lonband.WriteArray(lon_array)
    
    ds = gdal.Open("/home/yad2/sarpixeltracking/geometry_tif/lat.tif")
    lat_array = ds.GetRasterBand(1).ReadAsArray()
    ds=None
    
    latband = outdata.GetRasterBand(3)
    latband.WriteArray(lat_array)
    
    outdata = None

    ###
       
    print("\n - Generating amplitude VRT file")
    amp_vrt_filename = os.path.join(
        os.getcwd(), file_name_root + "_amp", file_name_root + "_amp.vrt"
    )
    print(amp_vrt_filename)
    
    lon_file = os.path.join(os.getcwd(), "geometry_tif/lon.tif")
    
    lat_file = os.path.join(os.getcwd(), "geometry_tif/lat.tif")

    os.system(f"gdal_translate -of VRT {amp_filename} {amp_vrt_filename}")

    et = xml.etree.ElementTree.parse(amp_vrt_filename)
    
    geogcs = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]]'


    geoloc_tag = xml.etree.ElementTree.SubElement(et.getroot(), "Metadata")
    geoloc_tag.attrib["domain"] = "GEOLOCATION"
    
    mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    mdi.attrib["key"] = "SRS"
    mdi.text = geogcs

    # mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    # mdi.attrib["key"] = "X_DATASET"
    # mdi.text = lon_file

    mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    mdi.attrib["key"] = "X_BAND"
    mdi.text = "2"

    # mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    # mdi.attrib["key"] = "Y_DATASET"
    # mdi.text = lat_file

    mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    mdi.attrib["key"] = "Y_BAND"
    mdi.text = "3"

    mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    mdi.attrib["key"] = "PIXEL_OFFSET"
    mdi.text = "0"

    mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    mdi.attrib["key"] = "LINE_OFFSET"
    mdi.text = "0"

    mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    mdi.attrib["key"] = "PIXEL_STEP"
    mdi.text = "1"

    mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    mdi.attrib["key"] = "LINE_STEP"
    mdi.text = "1"
    
    mdi = xml.etree.ElementTree.SubElement(geoloc_tag, "mdi")
    mdi.attrib["key"] = "GEOREFERENCING_CONVENTION"
    mdi.text = "TOP_LEFT_CORNER"
    
    et.write(amp_vrt_filename)

    ###
    
    os.system(f"rm -rf {file_name_root}_reproj")
    os.mkdir(file_name_root + "_reproj")
    
    
    print("\n - Generating reprojected amplitude file")
    reproj_amp_filename = os.path.join(
        os.getcwd(), file_name_root + "_reproj", file_name_root + "_reproj.tif"
    )
    print(reproj_amp_filename)

    os.system(f"gdalwarp -geoloc -t_srs EPSG:4326 -overwrite {amp_filename} {reproj_amp_filename}")
    
def prepare_geolocation_rasters():
    print("\n - Modifying geometry files")
    
    lon_file = os.path.join(os.getcwd(), "geometry/lon.rdr.full")
    
    lat_file = os.path.join(os.getcwd(), "geometry/lat.rdr.full")
    
    os.system(f"rm -rf {os.path.join(os.getcwd(), 'geometry_tif')}")
    os.mkdir(os.path.join(os.getcwd(), "geometry_tif"))
    
    lon_tif_file = os.path.join(os.getcwd(), "geometry_tif", "lon.tif")
    
    lat_tif_file = os.path.join(os.getcwd(), "geometry_tif", "lat.tif")
    
    generate_simple_GeoTIFF(lon_file, lon_tif_file, flip_ud=True)
    
    generate_simple_GeoTIFF(lat_file, lat_tif_file, flip_ud=True)
    
    lon_tif_vrt_file = os.path.join(os.getcwd(), "geometry_tif", "lon.vrt")
    
    lat_tif_vrt_file = os.path.join(os.getcwd(), "geometry_tif", "lat.vrt")
    
    os.system(f"gdal_translate -of VRT {lon_tif_file} {lon_tif_vrt_file}")
    
    os.system(f"gdal_translate -of VRT {lat_tif_file} {lat_tif_vrt_file}")
    
    geogcs = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]]'
    
    et = xml.etree.ElementTree.parse(lon_tif_vrt_file)
    
    SRS_tag = xml.etree.ElementTree.SubElement(et.getroot(), "SRS")
    SRS_tag.text = geogcs
    
    et.write(lon_tif_vrt_file)
    
    et = xml.etree.ElementTree.parse(lat_tif_vrt_file)
    
    SRS_tag = xml.etree.ElementTree.SubElement(et.getroot(), "SRS")
    SRS_tag.text = geogcs
    
    et.write(lat_tif_vrt_file)
    




### Generate netCDF files that include geolocation rasters

prepare_geolocation_rasters()

file_name_root = "reference"

if os.path.exists("reference_slc_crop"):
    fn1 = "reference_slc_crop/reference.slc"
else:
    fn1 = "reference_slc/reference.slc"

in_ds = gdal.Open(fn1, gdal.GA_ReadOnly)
array1 = np.flipud(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
in_ds = None

generate_geocoded_tif(array1, file_name_root)


folder2 = "secondary_nc"
