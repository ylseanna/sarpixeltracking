def runAutoRIFT():
    import os
    
    if os.path.exists("reference_slc_crop"):
         os.system(
            "scripts/testautoRIFT.py -m reference_slc_crop/reference.slc -s secondary_slc_crop/secondary.slc"
        )
    else:
         os.system(
            "scripts/testautoRIFT.py -m reference_slc/reference.slc -s secondary_slc/secondary.slc"
        )

   

