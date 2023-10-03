def copy_to_dest(logger, destination):
    import json
    import os
    import glob

    print("Destination: "+destination)

    if os.path.exists(destination) == False:
        print("*** creating destination folder ***")
        os.mkdir(destination)

    
    log = open('log.json')
    logdata = json.load(log)

    reference = logdata["frame_metadata"][0]["reference"]
    secondary = logdata["frame_metadata"][1]["secondary"]
    
    subfolder = reference['begintime'][:10].replace("-", "")+'-'+secondary['begintime'][:10].replace("-", "")

    print("\nSubfolder: "+subfolder)

    if os.path.exists(os.path.join(destination, subfolder)) == False:
        print("*** creating subfolder ***")
        os.mkdir(os.path.join(destination, subfolder))

    to_copy = [
        "preview",
        "geocoded_offsets",
        "log.json",
        "isce.log",
        "offset.mat"
    ]

    for xml in glob.glob("*.xml"):
        to_copy.append(xml)

    print("\nCopying files: ")

    for origin_folder in to_copy:
        if os.path.exists(origin_folder):
            os.system(f"cp -Rv {origin_folder} {os.path.join(destination, subfolder)}")
