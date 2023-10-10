def runISCE(logger, inps):
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

    os.system("stripmapApp.py stripmapApp.xml --start=startup  --end=topo")
    
    ### DENSEOFFSETS

    if inps.interferogram == True or inps.denseOffsets == True:
        logger.log("offset_resampling_start", "Starting offsets and resampling")

        os.system(
            "stripmapApp.py stripmapApp.xml --start=geo2rdr --end=refined_resample"
        )

        logger.log("offset_resampling_end", "Offsets and resampling finished")
    
    ### DENSEOFFSETS

    if inps.denseOffsets == True or inps.interferogram == True:
        logger.log("denseOffsets_start", "Starting dense offsets")


        os.system(
            "stripmapApp.py stripmapApp.xml --start=dense_offsets --end=dense_offsets"
        )

        logger.log("denseOffsets_end", "Dense offsets finished")
    else:
        logger.log("denseOffsets_skip", "Dense offsets skipped...")
        
        
    if inps.interferogram == True:
        logger.log("interferogram_start", "Starting interferogram processing")

        os.system(
            "stripmapApp.py stripmapApp.xml --start=rubber_sheet_range --end=endup"
        )

        logger.log("interferogram_end", "Interferogram processing finished")
    else:
        logger.log("interferogram_skip", "Interferogram processing skipped...")

    # os.system(
    #         "stripmapApp.py stripmapApp.xml --start=startup --end=formslc"
    # )

    # os.system(
    #         "stripmapApp.py stripmapApp.xml --start=verifyDEM --end=refined_resample"
    # )
