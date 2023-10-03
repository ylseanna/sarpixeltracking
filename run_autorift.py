def runAutoRIFT():
    import os

    os.system(
        "scripts/testautoRIFT.py -m reference_slc_crop/reference.slc -s coregisteredSlc/refined_coreg.slc"
    )

