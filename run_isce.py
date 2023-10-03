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
