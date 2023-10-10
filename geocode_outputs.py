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


def generate_previews():
    print(" - Geocoding and subsampling .slc images\n")
    from osgeo import gdal
    import numpy as np

    Geometry = generateGeometry()

    import os

    os.system("rm -rf preview")

    os.mkdir("preview")

    print("\nGeocoding and subsampling reference .slc image using GDAL:\n")
    
    
    
    if os.path.exists("reference_slc_crop"):
        in_filename = "reference_slc_crop/reference.slc"
    else:
        in_filename = "reference_slc/reference.slc"

    in_ds = gdal.Open(in_filename, gdal.GA_ReadOnly)
    in_array = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
    in_ds = None

    out_filename = "reference_preview"

    generateGeotiff(in_array, out_filename, "preview", Geometry, downsample=True)

    print("\nGeocoding and subsampling secondary .slc image using GDAL:\n")
    
    if os.path.exists("secondary_slc_crop"):
        in_filename = "secondary_slc_crop/secondary.slc"
    else:
        in_filename = "secondary_slc/secondary.slc"


    in_ds = gdal.Open(in_filename, gdal.GA_ReadOnly)
    in_array = np.fliplr(np.abs(in_ds.GetRasterBand(1).ReadAsArray()))
    in_ds = None

    out_filename = "secondary_preview"

    generateGeotiff(in_array, out_filename, "preview", Geometry, downsample=True)

    # ISN93 EPSG:3057


def geocode_offsets(inps):
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
            

def geocode_autoRIFT(inps):
    """This function specifically geocodes the .mat file and transfers it to an xyz type table for use in plotting, implicitly uses the ISCE geometry files"""
    from osgeo import gdal
    import numpy as np

    print(" - Importing .mat file\n")

    try:
        import h5py

        offset = h5py.File("offset.mat", "r")
    except:
        import scipy.io as sio

        offset = sio.loadmat("offset.mat")

    dx = offset["Dx"]
    dy = offset["Dy"]

    print(" - Reading ISCE geometry files\n")

    lats = gdal.Open("geometry/lat.rdr.full")
    lons = gdal.Open("geometry/lon.rdr.full")
    zs = gdal.Open("geometry/z.rdr.full")
    loss = gdal.Open("geometry/los.rdr.full")

    lat = np.array(lats.GetRasterBand(1).ReadAsArray())
    lon = np.array(lons.GetRasterBand(1).ReadAsArray())
    z = np.array(zs.GetRasterBand(1).ReadAsArray())
    los = np.array(loss.GetRasterBand(1).ReadAsArray())
    head = np.array(loss.GetRasterBand(2).ReadAsArray())

    lats = None
    lons = None
    zs = None
    loss = None

    off_height, off_width = np.shape(dx)

    isce_height, isce_width = np.shape(lat)

    print("AutoRIFT dimensions:", off_height, off_width)

    print("ISCE dimensions:", isce_height, isce_width)

    height_conv = isce_height / off_height
    width_conv = isce_width / off_width

    print("\n - Interpolating ISCE geometry files for application to AutoRIFT offsets\n")

    def interpolate(array, isce_width, isce_height, off_width, off_height):
        from scipy.interpolate import RegularGridInterpolator

        xrange = lambda x: np.linspace(0, 1, x)

        f = RegularGridInterpolator(
            (xrange(isce_width), xrange(isce_height)),
            array.T,
            method="linear",
            bounds_error=False,
        )

        xxnew, yynew = np.meshgrid(
            xrange(off_width), xrange(off_height), indexing="ij", sparse=True
        )

        interp = f((xxnew, yynew)).T

        return interp

    interp_lat = interpolate(lat, isce_width, isce_height, off_width, off_height)
    interp_lon = interpolate(lon, isce_width, isce_height, off_width, off_height)
    interp_z = interpolate(z, isce_width, isce_height, off_width, off_height)
    interp_los = interpolate(los, isce_width, isce_height, off_width, off_height)
    interp_head = interpolate(head, isce_width, isce_height, off_width, off_height)

    print(" - Saving data file\n")


    valid_vals = np.empty(shape=(0, 7))

    for i in range(off_height):
        for j in range(off_width):
            if not np.isnan(dx[i, j]):
                valid_vals = np.append(
                    valid_vals,
                    [
                        [
                            dx[i, j],
                            dy[i, j],
                            interp_lat[i, j],
                            interp_lon[i, j],
                            interp_z[i, j],
                            interp_los[i, j],
                            interp_head[i, j]
                        ]
                    ],
                    axis=0,
                )

    print(f"Number of valid autoRIFT values: {len(valid_vals)}")

    np.savetxt(
        "geocoded_offsets/AutoRIFT.data",
        valid_vals,
        fmt="%25.15f",
        delimiter=",",
        header="Dx, Dy, Lat, Lon, z, incidence, heading",
    )
