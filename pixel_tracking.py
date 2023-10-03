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
        help="Determines whether to ignore Ampcor standard offsets as generated for the ISCE program in further processing.",
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
        from geocode_outputs import geocode_offsets
        
        logger.log("geocode_start", "Starting geocoding offsets")

        geocode_offsets(inps)

        logger.log("geocode_end", "Geocoding offsets finished")
    else:
        logger.log("geocode_skip", "Geocoding offsets skipped...")

    ### Generate previews

    if inps.preview == True:
        from geocode_outputs import generate_previews
        
        logger.log("previews_start", "Starting generating previews")

        generate_previews()

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
