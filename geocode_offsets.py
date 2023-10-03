# SET UP IMPORTS


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
