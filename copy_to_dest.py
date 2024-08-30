def copy_to_dest(logger, destination, copyDestPre):
    import json
    import os
    from glob import glob

    print("Destination: "+destination)

    if os.path.exists(destination) == False:
        print("*** creating destination folder ***")
        os.mkdir(destination)

    
    log = open('log.json')
    logdata = json.load(log)

    reference = logdata["frame_metadata"][0]["reference"]
    secondary = logdata["frame_metadata"][1]["secondary"]
    
    if copyDestPre == None:
        subfolder = reference['begintime'][:10].replace("-", "")+'-'+secondary['begintime'][:10].replace("-", "")
    else:
        subfolder = copyDestPre+reference['begintime'][:10].replace("-", "")+'-'+secondary['begintime'][:10].replace("-", "")
    
    
    
    # subfolder = "test_folder"

    print("\nSubfolder: "+subfolder)

    if os.path.exists(os.path.join(destination, subfolder)) == False:
        print("*** creating subfolder ***")
        os.mkdir(os.path.join(destination, subfolder))
        
        os.mkdir(os.path.join(destination, subfolder, 'vmap'))
    
    

    to_copy_base = [
        "log.json",
        "isce.log",
        "vmap-offsets.csv"
    ]

    for xml in glob("*.xml"):
        to_copy_base.append(xml)

    for tif in glob("*_tif"):
        to_copy_base.append(tif)
        
    print(to_copy_base)

    vmap_dir = glob("reference*vmap*")[0]
    
    to_copy_vmap = [
        vmap_dir+"/vmap-F.tif",
        vmap_dir+"/vmap-dx.tif",
        vmap_dir+"/vmap-dy.tif"
    ]
    
    for tif in glob(vmap_dir+"/*.txt"):
        to_copy_vmap.append(tif)
        
    print(to_copy_vmap)

    print("\nCopying files: ")

    for file_or_folder in to_copy_base:
        if os.path.exists(file_or_folder):
            os.system(f"cp -Rv {file_or_folder} {os.path.join(destination, subfolder)}")


    for file_or_folder in to_copy_vmap:
        if os.path.exists(file_or_folder):
            os.system(f"cp -Rv {file_or_folder} {os.path.join(destination, subfolder, 'vmap')}")