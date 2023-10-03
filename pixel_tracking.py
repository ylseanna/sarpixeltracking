#!/usr/bin/env python

"""
Base function for pixel tracking for this project

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
            "log_metadata": {
                "program": "Pixel Tracking Pipeline",
                "program_desc": "Generates pixel tracking offset in SAR data for various platforms",
                "log_timestamp": t.timestamp(),
                "log_isotime": t.isoformat(),
            },
            "frame_metadata": [],
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
    
    def addFrameMetadata(self, frame_name, frame_data):
        import json

        logdata = {
            frame_name : frame_data,
        }

        with open("log.json", "r+") as file:
            json_object = json.load(file)

            json_object["frame_metadata"].append(logdata)
            # Sets file's current position at offset.
            file.seek(0)
            # convert back to json.
            json.dump(json_object, file, indent=4)


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


### Function for generating geotiffs based on geometry and array provided


def generateGeotiff(
    array, out_filename, folder, geometry, **kwargs
):  # no extension on filename
    downsample = kwargs.get(
        "downsample", False
    )  # optional arguments for after gdalwapr

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

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    out_ds.SetProjection(srs.ExportToWkt())

    band = out_ds.GetRasterBand(1)
    band.WriteArray(array)
    band.FlushCache()
    band.ComputeStatistics(False)

    out_ds = None

    os.system(
        f"gdal_translate -r bilinear\
                            -gcp 0 0 {geometry['UpperLeft'][1]} {geometry['UpperLeft'][0]} \
                            -gcp {width} 0 {geometry['UpperRight'][1]} {geometry['UpperRight'][0]} \
                            -gcp {width} {height} {geometry['LowerRight'][1]} {geometry['LowerRight'][0]} \
                            -gcp 0 {height} {geometry['LowerLeft'][1]} {geometry['LowerLeft'][0]} \
                            {os.path.join(folder, out_filename+'_untranslated.tif')} \
                            {os.path.join(folder, out_filename+'_unwarped.tif')}"
    )

    if downsample == True:
        out_name = os.path.join(folder, out_filename+'_uncompressed.tif')
    else:
        out_name = os.path.join(folder, out_filename+'.tif')

    os.system(
        f"gdalwarp \
                   -r bilinear -t_srs EPSG:4326 -et 0 -dstnodata nan \
                  {os.path.join(folder, out_filename+'_unwarped.tif')} \
                  {out_name}"
    )

    if downsample == True:
        print("\nExtra downsampling step")

        os.system(
            f"gdal_translate -r bilinear \
                -outsize 3840 0 \
                -co COMPRESS=LERC_ZSTD \
                {os.path.join(folder, out_filename+'_uncompressed.tif')} \
                {os.path.join(folder, out_filename+'.tif')}"
        )

        os.system

    os.system(f"rm -rf {os.path.join(folder, out_filename+'_untranslated.tif')}")
    os.system(f"rm -rf {os.path.join(folder, out_filename+'_unwarped.tif')}")
    if downsample == True:
        os.system(f"rm -rf {os.path.join(folder, out_filename+'_uncompressed.tif')}")


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
        required=True,
        help="Input filename 1 as found in the metadata database, typically reference. The program will automatically assign the earlier images as reference.",
    )
    parser.add_argument(
        "-f2",
        "--file2",
        dest="file2",
        type=str,
        required=True,
        help="Input filename 2 as found in the metadata database, typically secondary.",
    )
    parser.add_argument(
        "-d",
        "--dem",
        dest="DEM_file",
        type=str,
        required=True,
        help=".dem file to be used by ISCE for georeferencing",
    )
    parser.add_argument(
        "-dest",
        "--destination",
        dest="destination",
        type=str,
        required=False,  # True
        default=None,
        help="Target folder to store pipeline results (logs, previews, geocoded offsets), if this parameter is not set, no files will be copied",
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
        "--ignore-ampcor",
        default=False,
        action=argparse.BooleanOptionalAction,
        dest="ignore_ampcor",
        required=False,
        help="Determines whether to ignore Ampcor standard offsets as generated for the ISCE program in further processing, e.g.: geocoding (NOTE: they will always be generated).",
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
    parser.add_argument(
        "--previews",
        default=True,
        action=argparse.BooleanOptionalAction,
        dest="preview",
        required=False,
        help="Determines whether to generate previews of the .slc images for later visualisation (they are significantly downsampled)",
    )
    return parser.parse_args()


def generatePreviews():
    print(" - Geocoding and subsampling .slc images\n")
    from osgeo import gdal
    import numpy as np

    Geometry = generateGeometry()

    import os

    os.system("rm -rf preview")

    os.mkdir("preview")

    print("\nGeocoding and subsampling reference .slc image using GDAL:\n")

    in_filename = "reference_slc_crop/reference.slc"

    in_ds = gdal.Open(in_filename, gdal.GA_ReadOnly)
    in_array = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
    in_ds = None

    out_filename = "reference_preview"

    generateGeotiff(in_array, out_filename, "preview", Geometry, downsample=True)

    print("\nGeocoding and subsampling secondary .slc image using GDAL:\n")

    in_filename = "coregisteredSlc/refined_coreg.slc"

    in_ds = gdal.Open(in_filename, gdal.GA_ReadOnly)
    in_array = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
    in_ds = None

    out_filename = "secondary_preview"

    generateGeotiff(in_array, out_filename, "preview", Geometry, downsample=True)

    # ISN93 EPSG:3057


def geocodeOffsets(inps):
    import os
    import os.path
    from osgeo import gdal
    import numpy as np

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

    if inps.ignore_ampcor == True:
        print("\n*** Ignoring ampcor offsets (improves speed) ***")

    ### finding geometry

    Geometry = generateGeometry()

    # cleanup and intialize folder

    print("\n - Starting geocode\n")

    os.system("rm -rf geocoded_offsets")

    os.mkdir("geocoded_offsets")

    # geocoding

    for geocode in to_geocode:
        if geocode["program"] == "autoRIFT":
            os.system(f"rm -rf {os.path.join('geocoded_offsets', geocode['program'])}")
            os.mkdir(os.path.join("geocoded_offsets", geocode["program"]))

            print(f"Geocoding for {geocode['program']}\n")

            try:
                import h5py

                f = h5py.File("offset.mat", "r")
            except:
                import scipy.io as sio

                f = sio.loadmat("offset.mat")

            print(f"Geocoding pixel range offset:\n")

            xarray = np.fliplr(f["Dx"])

            generateGeotiff(
                xarray,
                "range_radar",
                os.path.join("geocoded_offsets", geocode["program"]),
                Geometry,
            )

            print(f"\nPixel to geographic distance conversion:")

            xconv = (
                (Geometry["top_geodesic"]["s12"] + Geometry["bottom_geodesic"]["s12"])
                / 2
                / Geometry["pixel_width"]
            )
            print(f"{xconv} meter/pixel\n")

            xarray_conv = xarray * xconv

            generateGeotiff(
                xarray_conv,
                "range",
                os.path.join("geocoded_offsets", geocode["program"]),
                Geometry,
            )

            print(f"\nGeocoding pixel azimuth offset:\n")

            yarray = np.fliplr(f["Dy"])

            generateGeotiff(
                yarray,
                "azimuth_radar",
                os.path.join("geocoded_offsets", geocode["program"]),
                Geometry,
            )

            print(f"\nPixel to geographic distance conversion:")

            yconv = (
                (Geometry["left_geodesic"]["s12"] + Geometry["right_geodesic"]["s12"])
                / 2
                / Geometry["pixel_height"]
            )
            print(f"{yconv} meter/pixel\n")

            yarray_conv = yarray * yconv

            generateGeotiff(
                yarray_conv,
                "azimuth",
                os.path.join("geocoded_offsets", geocode["program"]),
                Geometry,
            )

        elif geocode["program"] == "Ampcor" and inps.ignore_ampcor == False:
            os.system(f"rm -rf {os.path.join('geocoded_offsets', geocode['program'])}")
            os.mkdir(os.path.join("geocoded_offsets", geocode["program"]))

            print(f"Geocoding for {geocode['program']}\n")

            print(f"Geocoding range offset:\n")

            in_ds = gdal.Open("offsets/range.off", gdal.GA_ReadOnly)
            xarray = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
            in_ds = None

            generateGeotiff(
                xarray,
                "range",
                os.path.join("geocoded_offsets", geocode["program"]),
                Geometry,
            )

            print(f"\nGeocoding azimuth offset:\n")

            in_ds = gdal.Open("offsets/azimuth.off", gdal.GA_ReadOnly)
            yarray = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
            in_ds = None

            generateGeotiff(
                yarray,
                "azimuth",
                os.path.join("geocoded_offsets", geocode["program"]),
                Geometry,
            )
        elif geocode["program"] == "DenseAmpcor" and inps.ignore_ampcor == False:
            os.system(f"rm -rf {os.path.join('geocoded_offsets', geocode['program'])}")
            os.mkdir(os.path.join("geocoded_offsets", geocode["program"]))

            print(f"Geocoding for {geocode['program']}\n")

            # print(f"Geocoding range offset:\n")

            # xarray = np.fliplr(f['Dx'])

            # generateGeotiff(xarray, "range", geocode['program'], Geometry)

            # print(f"\nGeocoding azimuth offset:\n")

            # yarray = np.fliplr(f['Dy'])

            # generateGeotiff(yarray, "azimuth", geocode['program'], Geometry)



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
        from init_pipeline import init
        
        logger.log("init_start", "Starting initialisation")

        init(logger, inps.file1, inps.file2, inps.DEM_file)

        logger.log("init_end", "Initialisation finished")
    else:
        logger.log("init_skip", "Initalisation skipped...")

    ### ISCE

    if inps.isce == True:
        from run_isce import runISCE
        
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


    ### Start autoRIFT

    if inps.autoRIFT == True:
        from run_autorift import runAutoRIFT
        
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

        geocodeOffsets(inps)

        logger.log("geocode_end", "Geocoding offsets finished")
    else:
        logger.log("geocode_skip", "Geocoding offsets skipped...")

    ### Generate previews

    if inps.preview == True:
        logger.log("previews_start", "Starting generating previews")

        generatePreviews()

        logger.log("previews_end", "Generating previews finished")
    else:
        logger.log("previews_skip", "Generating previews skipped...")

    ### Generate previews

    if inps.destination != None:
        from copy_to_dest import copy_to_dest
        
        logger.log("copy_start", "Copying results to destination folder")

        copy_to_dest(logger, inps.destination)

        logger.log("copy_end", "Copying finished")
    else:
        logger.log("copy_skip", "Copying results to destination folder skipped...")

    ### FINISH

    logger.log("program_end", "Processing finished")
    logger.end()


if __name__ == "__main__":
    main()
