#!/usr/bin/env python

"""
Base function for pixel tracking for this project

Requires a file called MetaData.db containing the file locations of all files used in this project

requires: gdal, geographiclib
"""

### utils


def timestamp(string):
    from datetime import datetime

    return datetime.strptime(string, "%Y-%m-%dT%H:%M:%S")


class Logger:
    def __init__(self, inps):
        import json
        from datetime import datetime

        t = datetime.now()

        self.starttime = t

        args = {}

        for key, value in vars(inps).items():
            args[key] = value

        jsondata = {
            "metadata": {
                "program": "Pixel Tracking Pipeline",
                "program_desc": "Generates pixel tracking offset in SAR data for various platforms",
                "log_timestamp": t.timestamp(),
                "log_isotime": t.isoformat(),
            },
            "args": args,
            "logs": [],
            "errors": [],
        }

        json_object = json.dumps(jsondata, indent=4)

        with open("log.json", "w") as outfile:
            outfile.write(json_object)

    def log(self, id, description):
        import json
        from datetime import datetime

        t = datetime.now()

        logdata = {
            "id": id,
            "level": "log",
            "timestamp": t.timestamp(),
            "isotime": t.isoformat(),
            "descr": description,
        }

        with open("log.json", "r+") as file:
            json_object = json.load(file)

            json_object["logs"].append(logdata)
            # Sets file's current position at offset.
            file.seek(0)
            # convert back to json.
            json.dump(json_object, file, indent=4)

        print(f"\n{t} - {description}\n")

    def error(self, exctype, value, traceback):
        import json
        from datetime import datetime

        t = datetime.now()

        errordata = {
            "exctype": exctype,
            "value": value,
            "traceback": traceback,
            "timestamp": t.timestamp(),
            "isotime": t.isoformat(),
        }

        with open("log.json", "r+") as file:
            json_object = json.load(file)

            json_object["errors"].append(errordata)
            # Sets file's current position at offset.
            file.seek(0)
            # convert back to json.
            json.dump(json_object, file, indent=4)

    def end(self):
        from datetime import datetime

        program_end = datetime.now()

        print(f"Completed in {program_end - self.starttime}")


def ErrorHandling(logger):
    def log_exception(exctype, value, traceback):
        logger.error(str(exctype), str(value), str(traceback))

    def attach_hook(hook_func, run_func):
        def inner(*args, **kwargs):
            if not (args or kwargs):
                # This condition is for sys.exc_info
                local_args = run_func()
                hook_func(*local_args)
            else:
                # This condition is for sys.excepthook
                hook_func(*args, **kwargs)
            return run_func(*args, **kwargs)

        return inner

    import sys

    sys.exc_info = attach_hook(log_exception, sys.exc_info)
    sys.excepthook = attach_hook(log_exception, sys.excepthook)


### Main functions


def argparse():
    """
    Command line parser. Adapted from testGeogrid.py
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="ERS Pixel Tracking pipeline",
        description="Completes a full run of the ERS pixel tracking pipelines for single frame (stripmap) images.",
    )
    parser.add_argument(
        "-f1",
        "--file1",
        dest="file1",
        type=str,
        required=False,  # Set to true eventually
        help="Input filename 1 as found in the metadata database, typically reference. The program will automatically assign the earlier images as reference.",
    )
    parser.add_argument(
        "-f2",
        "--file2",
        dest="file2",
        type=str,
        required=False,  # Set to true eventually
        help="Input filename 2 as found in the metadata database, typically secondary.",
    )
    parser.add_argument(
        "-d",
        "--dem",
        dest="DEM_file",
        type=str,
        required=False,  # True
        help="Input dem file as found in the metadata database.",
    )
    parser.add_argument(
        "--init",
        default=True,
        action=argparse.BooleanOptionalAction,
        dest="init",
        required=False,
        help="Determines whether to do go through the initialisation phase of the program, e.g. setting up files and folder structure etc. ",
    )
    parser.add_argument(
        "--isce",
        default=True,
        action=argparse.BooleanOptionalAction,
        dest="isce",
        required=False,
        help="Determines whether to run ISCE preprocessing step.",
    )
    parser.add_argument(
        "--denseOffsets",
        default=False,
        action=argparse.BooleanOptionalAction,
        dest="denseOffsets",
        required=False,
        help="Determines whether to do dense Ampcor dense offsets as provided in the ISCE program.",
    )
    parser.add_argument(
        "--Geogrid",
        default=False,
        action=argparse.BooleanOptionalAction,
        dest="geogrid",
        required=False,
        help="Determines whether to run Geogrid step (currently not working).",
    )
    parser.add_argument(
        "--autoRIFT",
        default=True,
        action=argparse.BooleanOptionalAction,
        dest="autoRIFT",
        required=False,
        help="Determines whether to run autoRIFT step.",
    )
    parser.add_argument(
        "--geocode",
        default=True,
        action=argparse.BooleanOptionalAction,
        dest="geocode",
        required=False,
        help="Determines whether to geocode the produced offset from each generator (ampcor, denseAmpcor, autoRIFT).",
    )
    return parser.parse_args()


def generatePreviews():
    print(" - Geocoding offsets\n")
    from osgeo import gdal
    import numpy as np

    print("Getting geometry bounds:")

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

    print("Upper left corner:", UpperLeft)
    print("Upper right corner:", UpperRight)
    print("Lower right corner", LowerRight)
    print("Lower left corner", LowerLeft)

    import os

    ### map project reference slc

    print("\nMap projecting the reference .slc image using GDAL\n")

    # cleanup

    os.system("rm -rf reference_geotiff")

    os.mkdir("reference_geotiff")

    out_filename = "reference_geotiff/reference.tif"

    in_filename = "reference_slc/reference.slc"

    driver = gdal.GetDriverByName("GTiff")

    in_ds = gdal.Open(in_filename, gdal.GA_ReadOnly)
    arr = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))

    out_ds = driver.Create(
        out_filename, arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32
    )
    out_ds.SetProjection(in_ds.GetProjection())

    band = out_ds.GetRasterBand(1)
    band.WriteArray(arr)
    band.FlushCache()
    band.ComputeStatistics(False)

    in_ds = None
    out_ds = None

    command = f"gdal_translate \
                    -gcp 0 0 {UpperLeft[1]} {UpperLeft[0]} \
                    -gcp {width} 0 {UpperRight[1]} {UpperRight[0]} \
                    -gcp {width} {height} {LowerRight[1]} {LowerRight[0]} \
                    -gcp 0 {height} {LowerLeft[1]} {LowerLeft[0]} \
                    reference_geotiff/reference.tif reference_geotiff/reference_projec.tif"

    os.system(command)

    os.system(
        "gdalwarp -r bilinear -t_srs EPSG:32627 -et 0 -dstnodata -32767 reference_geotiff/reference_projec.tif reference_geotiff/reference_warp.tif"
    )

    ### map project secondary slc

    print("\nMap projecting the coregistered secondary .slc image using GDAL\n")

    # cleanup

    os.system("rm -rf secondary_geotiff")

    os.mkdir("secondary_geotiff")

    out_filename = "secondary_geotiff/secondary.tif"

    in_filename = "coregisteredSlc/refined_coreg.slc"

    driver = gdal.GetDriverByName("GTiff")

    in_ds = gdal.Open(in_filename, gdal.GA_ReadOnly)
    arr = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))

    out_ds = driver.Create(
        out_filename, arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32
    )
    out_ds.SetProjection(in_ds.GetProjection())

    band = out_ds.GetRasterBand(1)
    band.WriteArray(arr)
    band.FlushCache()
    band.ComputeStatistics(False)

    in_ds = None
    out_ds = None

    command = f"gdal_translate \
                    -gcp 0 0 {UpperLeft[1]} {UpperLeft[0]} \
                    -gcp {width} 0 {UpperRight[1]} {UpperRight[0]} \
                    -gcp {width} {height} {LowerRight[1]} {LowerRight[0]} \
                    -gcp 0 {height} {LowerLeft[1]} {LowerLeft[0]} \
                    secondary_geotiff/secondary.tif secondary_geotiff/secondary_projec.tif"

    os.system(command)

    os.system(
        "gdalwarp -r bilinear -t_srs EPSG:32627 -et 0 -dstnodata -32767 secondary_geotiff/secondary_projec.tif secondary_geotiff/secondary_warp.tif"
    )

    # ISN93 EPSG:3057


def geocodeOffsets():
    import os
    import os.path
    from osgeo import gdal, osr
    import numpy as np

    ### Function for generating geotiffs based on geometry and array provided

    def generateGeotiff(array, out_filename, folder, geometry):  # no extension on filename
        driver = gdal.GetDriverByName("GTiff")

        out_ds = driver.Create(
            os.path.join("geocoded_offsets", folder, out_filename+"_untranslated.tif"),
            array.shape[1],
            array.shape[0],
            1,
            gdal.GDT_Float32,
        )
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        out_ds.SetProjection(srs.ExportToWkt())

        band = out_ds.GetRasterBand(1)
        band.WriteArray(array)
        band.FlushCache()
        band.ComputeStatistics(False)

        out_ds = None

        os.system(f"gdal_translate \
                        -gcp 0 0 {geometry['UpperLeft'][1]} {geometry['UpperLeft'][0]} \
                        -gcp {width} 0 {geometry['UpperRight'][1]} {geometry['UpperRight'][0]} \
                        -gcp {width} {height} {geometry['LowerRight'][1]} {geometry['LowerRight'][0]} \
                        -gcp 0 {height} {geometry['LowerLeft'][1]} {geometry['LowerLeft'][0]} \
                        {os.path.join('geocoded_offsets', folder, out_filename+'_untranslated.tif')} \
                        {os.path.join('geocoded_offsets', folder, out_filename+'_unwarped.tif')}")

        os.system(f"gdalwarp \
                   -r bilinear -t_srs EPSG:3057 -et 0 \
                  {os.path.join('geocoded_offsets', folder, out_filename+'_unwarped.tif')} \
                  {os.path.join('geocoded_offsets', folder, out_filename+'.tif')}")

    ### setting up loop

    print(" - Detecting offsets to geocode\n")

    programs = [
        {"filepath": "offsets", "pathtype": "folder", "program": "Ampcor"},
        {"filepath": "denseOffsets", "pathtype": "folder", "program": "DenseAmpcor"},
        {"filepath": "offset.mat", "pathtype": ".mat file", "program": "autoRIFT"},
    ]

    to_geocode = []

    for program in programs:
        if os.path.exists(program["filepath"]):
            to_geocode.append(program)

            print(f"A {program['pathtype']} detected for {program['program']}.")

    ### finding geometry

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

    print("Upper left corner:", UpperLeft)
    print("Upper right corner:", UpperRight)
    print("Lower right corner", LowerRight)
    print("Lower left corner", LowerLeft)

    Geometry = {
        "UpperLeft": UpperLeft,
        "UpperRight": UpperRight,
        "LowerRight": LowerRight,
        "LowerLeft": LowerLeft,
    }

    # cleanup and intialize folder

    os.system("rm -rf geocoded_offsets")

    os.mkdir("geocoded_offsets")

    for geocode in to_geocode:
        os.system(f"rm -rf {os.path.join('geocoded_offsets', geocode['program'])}")
        os.mkdir(os.path.join("geocoded_offsets", geocode['program']))

        if geocode['program'] == 'autoRIFT':
            try:
                import h5py
                f = h5py.File('offset.mat','r')
            except:
                import scipy.io as sio
                f = sio.loadmat('offset.mat')

            array = np.fliplr(f['Dx'])

            print(array)

            generateGeotiff(array, "dx", geocode['program'], Geometry)

 



def runAutoRIFT():
    import os

    os.system(
        "scripts/testautoRIFT.py -m reference_slc/reference.slc -s coregisteredSlc/refined_coreg.slc"
    )


def runGeogridCompat(referenceFile, demFile):
    """
    Wire and run geogrid, taken from testGeogridOptical.py, but made to be compatible for the general translation step
    """

    from osgeo import gdal

    dem_info = gdal.Info(demFile, format="json")
    reference_info = gdal.Info(referenceFile, format="json")

    from geogrid import GeogridOptical

    obj = GeogridOptical()

    obj.startingX = reference_info["geoTransform"][0]
    obj.startingY = reference_info["geoTransform"][3]
    obj.XSize = reference_info["geoTransform"][1]
    obj.YSize = reference_info["geoTransform"][5]

    import xml.etree.ElementTree as ET
    from datetime import datetime

    time1 = (
        ET.parse("reference_slc.xml")
        .find(
            'component[@name="instance"]/component[@name="orbit"]/property[@name="min_time"]/value'
        )
        .text
    )
    time2 = (
        ET.parse("secondary_slc.xml")
        .find(
            'component[@name="instance"]/component[@name="orbit"]/property[@name="min_time"]/value'
        )
        .text
    )

    t1 = datetime.strptime(time1, "%Y-%m-%d %H:%M:%S.%f")
    t2 = datetime.strptime(time2, "%Y-%m-%d %H:%M:%S.%f")
    deltat = t2 - t1
    obj.repeatTime = deltat.total_seconds()

    obj.numberOfLines = reference_info["size"][1]
    obj.numberOfSamples = reference_info["size"][0]

    obj.nodata_out = -32767
    obj.chipSizeX0 = 240

    obj.gridSpacingX = dem_info["geoTransform"][1]

    obj.dat1name = referenceFile
    obj.demname = demFile
    #     obj.dhdxname = dhdx
    #     obj.dhdyname = dhdy
    #     obj.vxname = vx
    #     obj.vyname = vy
    #     obj.srxname = srx
    #     obj.sryname = sry
    #     obj.csminxname = csminx
    #     obj.csminyname = csminy
    #     obj.csmaxxname = csmaxx
    #     obj.csmaxyname = csmaxy
    # obj.ssmname = ssm
    obj.winlocname = "window_location.tif"
    obj.winoffname = "window_offset.tif"
    obj.winsrname = "window_search_range.tif"
    obj.wincsminname = "window_chip_size_min.tif"
    obj.wincsmaxname = "window_chip_size_max.tif"
    obj.winssmname = "window_stable_surface_mask.tif"
    obj.winro2vxname = "window_rdr_off2vel_x_vec.tif"
    obj.winro2vyname = "window_rdr_off2vel_y_vec.tif"
    obj.winsfname = "window_scale_factor.tif"
    ##dt-varying search range scale (srs) rountine parameters
    obj.srs_dt_unity = 32
    obj.srs_max_scale = 10
    obj.srs_max_search = 20000
    obj.srs_min_search = 0

    obj.runGeogrid()

    # run_info = {
    #     'chipsizex0': obj.chipSizeX0,
    #     'gridspacingx': obj.gridSpacingX,
    #     'vxname': vx,
    #     'vyname': vy,
    #     'sxname': kwargs.get('dhdxs'),
    #     'syname': kwargs.get('dhdys'),
    #     'maskname': kwargs.get('sp'),
    #     'xoff': obj.pOff,
    #     'yoff': obj.lOff,
    #     'xcount': obj.pCount,
    #     'ycount': obj.lCount,
    #     'dt': obj.repeatTime,
    #     'epsg': kwargs.get('epsg'),
    #     'XPixelSize': obj.X_res,
    #     'YPixelSize': obj.Y_res,
    #     'cen_lat': obj.cen_lat,
    #     'cen_lon': obj.cen_lon,
    # }

    # return run_info


def runISCE():
    import os
    from os import path

    # cleanup

    folders = [
        "PICKLE",
        "misreg",
        "geometry",
        "coregisteredSlc",
        "offsets",
        "denseOffsets",
        "interferogram",
    ]

    for folder in folders:
        if path.isdir(folder):
            os.system(f"rm -rf {folder}")

    if path.exists("demLat*"):
        os.system("rm demLat*")

    os.system("stripmapApp.py stripmapApp.xml --start=startup --end=refined_resample")

    # os.system(
    #         "stripmapApp.py stripmapApp.xml --start=startup --end=formslc"
    # )

    # os.system(
    #         "stripmapApp.py stripmapApp.xml --start=verifyDEM --end=refined_resample"
    # )


def init(file1, file2):
    """function that prepares the pipeline"""

    ### Imports

    import os, glob
    import sqlite3
    from datetime import datetime

    ### INPUT

    # cleanup

    os.system("rm -rf reference*")
    os.system("rm -rf secondary*")
    os.system("rm *.xml")

    # retrieve arguments

    """ eventually change to input """

    file1 = "ER01_SAR_IM__0P_19920524T123500_19920524T123517_UPA_04477_0000.CEOS.tar.gz"
    file2 = "ER01_SAR_IM__0P_19930718T123508_19930718T123525_UPA_10489_0000.CEOS.tar.gz"

    print("Files selected:\n" + file1 + "\n" + file2)

    ### QUERY FILES and ASSIGN REF AND SEC

    print("\n - Querying files...")

    def query(query_str):
        conn = sqlite3.connect("MetaData.db")

        cursor = conn.execute(query_str)

        results = {
            "results": [
                dict(zip([column[0] for column in cursor.description], row))
                for row in cursor.fetchall()
            ]
        }

        conn.close()

        return results

    files = query(
        "SELECT * from frames WHERE filename = '{}' OR filename = '{}'".format(
            file1, file2
        )
    )

    file1, file2 = files["results"]

    print("\nFiles succesfully queried:\n")

    timedelta = timestamp(file1["begintime"]) - timestamp(file2["begintime"])

    if timedelta.total_seconds() < 0:
        reference = file1
        secondary = file2
    else:
        reference = file2
        secondary = file1

    print(
        "Reference:\n  Name:          "
        + reference["filename"]
        + "\n  Begin time:    "
        + reference["begintime"]
        + "\n  End time:      "
        + reference["endtime"]
        + "\n  Sensor:        "
        + reference["sensor"]
        + "\n  File location: "
        + reference["fileloc"]
    )
    print(
        "\nSecondary:\n  Name:          "
        + secondary["filename"]
        + "\n  Begin time:    "
        + secondary["begintime"]
        + "\n  End time:      "
        + secondary["endtime"]
        + "\n  Sensor:        "
        + secondary["sensor"]
        + "\n  File location: "
        + secondary["fileloc"]
    )

    timedelta = timestamp(secondary["begintime"]) - timestamp(reference["begintime"])

    print(
        "\nPair information:\n  Baseline:      "
        + str(timedelta)
        + "\n  Track:         "
        + str(reference["track"])
        + "\n  Frame:         "
        + str(reference["frame"])
        + ""
    )

    ### CREATE XML files and processing folders

    ### /home/data/orbits/ODR/ERS1
    ### /home/data/orbits/ODR/ERS2

    # Folders:

    print("\n - Generating folder structure...")

    # make sure ref_dir exists and is empty
    if os.path.exists("./reference") == False:
        os.mkdir("reference")
    else:
        files = glob.glob("./reference/*")
        for f in files:
            os.remove(f)

    # make sure sec_dir exists and is empty
    if os.path.exists("./secondary") == False:
        os.mkdir("secondary")
    else:
        files = glob.glob("./secondary/*")
        for f in files:
            os.remove(f)

    # Unpack images:

    print("\n - Unpacking images...")

    print("\nReference:\n")

    os.system(f"tar -zvxf {reference['fileloc']} --directory ./reference")

    print("\nSecondary:\n")

    os.system(f"tar -zvxf {secondary['fileloc']} --directory ./secondary")

    # Generate XML-files

    print("\n - Generating XML-files...")

    # reference
    if reference["sensor"] == "ERS 1":
        reference_orbitloc = "/home/data/orbits/ODR/ERS1"
    elif reference["sensor"] == "ERS 2":
        reference_orbitloc = "/home/data/orbits/ODR/ERS2"

    print("\nreference.xml:")

    reference_xml = f"""
    <component name="Reference">
        <property name="IMAGEFILE">
            ./reference/DAT_01.001
        </property>
        <property name="LEADERFILE">
            ./reference/LEA_01.001
        </property>
        <property name="OUTPUT">reference</property>
        <property name="ORBIT_TYPE">
            <value>ODR</value>
        </property>
        <property name="ORBIT_DIRECTORY">
            <value>{reference_orbitloc}</value>
        </property>
    </component>"""

    print(reference_xml)

    f = open("reference.xml", "w")
    f.write(reference_xml)
    f.close()

    # secondary
    if secondary["sensor"] == "ERS 1":
        secondary_orbitloc = "/home/data/orbits/ODR/ERS1"
    elif secondary["sensor"] == "ERS 2":
        secondary_orbitloc = "/home/data/orbits/ODR/ERS2"

    print("\nsecondary.xml:")

    secondary_xml = f"""
    <component name="Secondary">
        <property name="IMAGEFILE">
            ./secondary/DAT_01.001
        </property>
        <property name="LEADERFILE">
            ./secondary/LEA_01.001
        </property>
        <property name="OUTPUT">secondary</property>
        <property name="ORBIT_TYPE">
            <value>ODR</value>
        </property>
        <property name="ORBIT_DIRECTORY">
            <value>{secondary_orbitloc}</value>
        </property>
    </component>"""

    print(secondary_xml)

    f = open("secondary.xml", "w")
    f.write(secondary_xml)
    f.close()

    # stripmapApp.xml
    print("\nstripmapApp.xml:")

    DEM_loc = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem"  # "/home/yad2/DEM/IslandsDEMv1.0_2x2m_zmasl_isn93_SouthMerge.tif"

    stripmapApp_xml = f"""
    <stripmapApp>
        <component name="insar">
            <property  name="Sensor name">ERS</property>
            <component name="reference">
                <catalog>reference.xml</catalog>
            </component>
            <component name="secondary">
                <catalog>secondary.xml</catalog>
            </component>
            <property name="demFilename">{DEM_loc}</property>
            <property name="do denseoffsets">True</property>
        </component>
    </stripmapApp>"""

    print(stripmapApp_xml)

    f = open("stripmapApp.xml", "w")
    f.write(stripmapApp_xml)
    f.close()


### Main loop


def main():
    """
    Main driver.
    """

    from datetime import datetime
    import os

    print("\n*** ESR 1&2 - Pixel Tracking Pipeline ***\n")

    ### COMMAND LINE ARGUMENTS

    inps = argparse()

    logger = Logger(inps)

    ErrorHandling(logger)

    # mark starttime

    logger.log("start", "Starting program")

    ### INITIALISATION

    if inps.init == True:
        logger.log("init_start", "Starting initialisation")

        init(inps.file1, inps.file2)

        logger.log("init_end", "Initialisation finished")
    else:
        logger.log("init_skip", "Initalisation skipped...")

    ### ISCE

    if inps.isce == True:
        logger.log("isce_start", "Starting ISCE")

        runISCE()

        logger.log("isce_end", "ISCE finished")
    else:
        logger.log("ISCE_skip", "ISCE skipped...")

    ### DENSEOFFSETS

    if inps.denseOffsets == True:
        logger.log("denseOffsets_start", "Starting dense offsets")

        os.system(
            "stripmapApp.py stripmapApp.xml --start=dense_offsets --end=dense_offsets"
        )

        logger.log("isce_end", "Dense offsets finished")

    ### Start Geogrid

    if inps.geogrid == True:
        logger.log("Geogrid_start", "Starting Geogrid")

        ## gdalwarp -r bilinear -t_srs EPSG:3057 -of GTiff demLat_N63_N65_Lon_W022_W018.dem.wgs84 demLat_N63_N65_Lon_W022_W018.dem.ISN93
        ## EPSG:32627

        runGeogridCompat(
            "reference_geotiff/reference_warp.tif",
            "/home/yad2/DEM/IslandsDEMv1.0_2x2m_zmasl_isn93_SouthMerge.tif",
        )

        # os.system(
        #     "testGeogridOptical.py -m reference_geotiff/reference_warp.tif -s secondary_geotiff/secondary_warp.tif -d demLat_N63_N65_Lon_W022_W018.dem.ISN93"
        # )

    ### Start autoRIFT

    if inps.autoRIFT == True:
        logger.log("autoRIFT_start", "Starting AutoRIFT")

        ## gdalwarp -r bilinear -t_srs EPSG:3057 -of GTiff demLat_N63_N65_Lon_W022_W018.dem.wgs84 demLat_N63_N65_Lon_W022_W018.dem.ISN93
        ## EPSG:32627

        runAutoRIFT()

        # os.system(
        #     "testGeogridOptical.py -m reference_geotiff/reference_warp.tif -s secondary_geotiff/secondary_warp.tif -d demLat_N63_N65_Lon_W022_W018.dem.ISN93"
        # )

        logger.log("autoRIFT_end", "AutoRIFT finished")
    else:
        logger.log("autoRIFT_skip", "AutoRIFT skipped...")

    ### Geocode offsets

    if inps.geocode == True:
        logger.log("geocode_start", "Starting geocoding offsets")

        geocodeOffsets()

        logger.log("geocode_end", "Geocoding offsets finished")
    else:
        logger.log("geocode_skip", "Geocoding offsets skipped...")

    ### FINISH

    logger.log("program_end", "Processing finished")
    logger.end()


if __name__ == "__main__":
    main()
