def timestamp(string):
    from datetime import datetime

    return datetime.strptime(string, "%Y-%m-%dT%H:%M:%S")


def init(logger, file1, file2, demfile):
    """function that prepares the pipeline"""

    ### Imports

    import os, glob
    from datetime import datetime

    ### INPUT

    # cleanup

    os.system("rm -rf reference*")
    os.system("rm -rf secondary*")
    os.system("rm -rf coreg_secondary*")
    os.system("rm -rf geometry*")
    os.system("rm -rf preview")
    os.system("rm *.xml")

    print("Files selected:\n" + file1 + "\n" + file2)
    
    ### QUERY FILES and ASSIGN REF AND SEC

    print("\n - Querying files...")

    if file1.endswith(".CEOS.tar.gz"):
        print("\nERS 1,2 or EnVISAT files in old raw data format detected\n")
        
        generate_files_ERS_CEOS(logger, file1, file2, demfile) # NEEDS TO BE ADJUSTED TO NOT USE THE DATABASE
    if file1.endswith((".E1", ".E2", ".N1")):
        print("\nEnVISAT data format detected\n")

        generate_files_Envisat_format(logger, file1, file2, demfile)
    if os.path.basename(file1).startswith(("TSX", "TDX")):
        print("TerraSAR-X data format detected\n")
        
        generate_files_TSX_format(logger, file1, file2, demfile)

def generate_files_ERS_CEOS(logger, file1, file2, demfile):
    import os
    from datetime import datetime
    from glob import glob
    

    reference = {}
    secondary = {}

    reference["fileloc"] = file1
    secondary["fileloc"] = file2

    reference["filename"] = os.path.basename(file1)
    secondary["filename"] = os.path.basename(file2)

    basef1, extf1 = os.path.splitext(os.path.basename(file1))

    basef2, extf2 = os.path.splitext(os.path.basename(file2))

    if basef1.startswith("ER01"):
        reference["platform"] = "ERS-1"
    elif basef1.startswith("ER02"):
        reference["platform"] = "ERS-2"

    if basef2.startswith("ER01"):
        secondary["platform"] = "ERS-1"
    elif basef2.startswith("ER02"):
        secondary["platform"] = "ERS-2"
        
    from pygeotools.lib import timelib


    reference["begintime"] = timelib.fn_getdatetime_list(basef1)[0]

    secondary["begintime"] = timelib.fn_getdatetime_list(basef2)[0]
    
    print(
        "Reference:\n  Name:          "
        + reference["filename"]
        + "\n  Begin time:    "
        + reference["begintime"].isoformat()
        + "\n  Platform:      "
        + reference["platform"]
        + "\n  File location: "
        + reference["fileloc"]
    )
    print(
        "\nSecondary:\n  Name:          "
        + secondary["filename"]
        + "\n  Begin time:    "
        + secondary["begintime"].isoformat()
        + "\n  Platform:      "
        + secondary["platform"]
        + "\n  File location: "
        + secondary["fileloc"]
    )

    timedelta = secondary["begintime"] - reference["begintime"]

    print(
        "\nPair information:\n  Baseline:      "
        + str(timedelta)
    )

    reference["begintime"] = reference["begintime"].isoformat()
    secondary["begintime"] = secondary["begintime"].isoformat()

    logger.addFrameMetadata("reference", reference)
    logger.addFrameMetadata("secondary", secondary)
    
    # Folders:

    print("\n - Generating folder structure...")

    # make sure ref_dir exists and is empty
    if os.path.exists("./reference") == False:
        os.mkdir("reference")
    else:
        files = glob("./reference/*")
        for f in files:
            os.remove(f)

    # make sure sec_dir exists and is empty
    if os.path.exists("./secondary") == False:
        os.mkdir("secondary")
    else:
        files = glob("./secondary/*")
        for f in files:
            os.remove(f)
            
            
    print("\n - Unpacking images...")

    print("\nReference:\n")

    os.system(f"tar -zvxf {reference['fileloc']} --directory ./reference")

    print("\nSecondary:\n")

    os.system(f"tar -zvxf {secondary['fileloc']} --directory ./secondary")
    
    

    # Generate XML-files

    print("\n - Generating XML-files...")

    # reference
    if reference["platform"] == "ERS-1":
        reference_orbitloc = "/home/data/orbits/ODR/ERS1"
    elif reference["platform"] == "ERS-2":
        reference_orbitloc = "/home/data/orbits/ODR/ERS2"

    print("\nreference.xml:")

    reference_xml = f"""
    <component name="Reference">
        <property name="IMAGEFILE">
            ./reference/DAT_01.001
        </property>
        <property name="LEADERFILE">
            ./reference/LEA_01.001
        </property>
        <property name="OUTPUT">reference</property>
        <property name="ORBIT_TYPE">
            <value>ODR</value>
        </property>
        <property name="ORBIT_DIRECTORY">
            <value>{reference_orbitloc}</value>
        </property>
    </component>"""

    print(reference_xml)

    f = open("reference.xml", "w")
    f.write(reference_xml)
    f.close()

    # secondary
    if secondary["platform"] == "ERS-1":
        secondary_orbitloc = "/home/data/orbits/ODR/ERS1"
    elif secondary["platform"] == "ERS-2":
        secondary_orbitloc = "/home/data/orbits/ODR/ERS2"

    print("\nsecondary.xml:")

    secondary_xml = f"""
    <component name="Secondary">
        <property name="IMAGEFILE">
            ./secondary/DAT_01.001
        </property>
        <property name="LEADERFILE">
            ./secondary/LEA_01.001
        </property>
        <property name="OUTPUT">secondary</property>
        <property name="ORBIT_TYPE">
            <value>ODR</value>
        </property>
        <property name="ORBIT_DIRECTORY">
            <value>{secondary_orbitloc}</value>
        </property>
    </component>"""

    print(secondary_xml)

    f = open("secondary.xml", "w")
    f.write(secondary_xml)
    f.close()

    # stripmapApp.xml
    print("\nstripmapApp.xml:")

    # DEM_loc = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem"  # "/home/yad2/DEM/IslandsDEMv1.0_2x2m_zmasl_isn93_SouthMerge.tif"

    stripmapApp_xml = f"""
    <stripmapApp>
        <component name="insar">
            <property  name="Sensor name">ERS</property>
            <component name="reference">
                <catalog>reference.xml</catalog>
            </component>
            <component name="secondary">
                <catalog>secondary.xml</catalog>
            </component>
            <property name="demFilename">{demfile}</property>
            <property name="do denseoffsets">True</property>
            <property name="regionOfInterest">[63.699855,63.583704,-19.476357,-19.205132]</property>
        </component>
    </stripmapApp>"""

    print(stripmapApp_xml)

    f = open("stripmapApp.xml", "w")
    f.write(stripmapApp_xml)
    f.close()

    # stripmapApp.xml
    print("\ndense.xml:")

    # DEM_loc = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem"  # "/home/yad2/DEM/IslandsDEMv1.0_2x2m_zmasl_isn93_SouthMerge.tif"

    dense_xml = f"""
    <stripmapAppDenseAmpcor>
        <component name="dense">
            <property name="Ampcor window width">64</property>
            <property name="Ampcor window height">256</property>
            <property name="Ampcor search window width">10</property>
            <property name="Ampcor search window height">40</property>
            <property name="Ampcor skip width">128</property>
            <property name="Ampcor skip height">32</property>
        </component>
    </stripmapAppDenseAmpcor>"""

    print(dense_xml)

    f = open("dense.xml", "w")
    f.write(dense_xml)
    f.close()

def generate_files_Envisat_format(logger, file1, file2, demfile):
    import os
    from datetime import datetime

    reference = {}
    secondary = {}

    reference["fileloc"] = file1
    secondary["fileloc"] = file2

    reference["filename"] = os.path.basename(file1)
    secondary["filename"] = os.path.basename(file2)

    basef1, extf1 = os.path.splitext(os.path.basename(file1))

    basef2, extf2 = os.path.splitext(os.path.basename(file2))

    if extf1 == ".E1":
        reference["platform"] = "ERS-1"
    elif extf1 == ".E2":
        reference["platform"] = "ERS-2"
    elif extf1 == ".N1":
        reference["platform"] = "Envisat"

    if extf2 == ".E1":
        secondary["platform"] = "ERS-1"
    elif extf2 == ".E2":
        secondary["platform"] = "ERS-2"
    elif extf2 == ".N1":
        secondary["platform"] = "Envisat"

    reference["track"], secondary["track"] = (
        int(basef1.split("_")[6]),
        int(basef2.split("_")[6]),
    )

    reference["orbit"], secondary["orbit"] = (
        int(basef1.split("_")[7]),
        int(basef2.split("_")[7]),
    )

    reference["begintime"] = datetime(
        year=int(basef1.split("_")[3][6:10]),
        month=int(basef1.split("_")[3][10:12]),
        day=int(basef1.split("_")[3][12:14]),
        hour=int(basef1.split("_")[4][0:2]),
        minute=int(basef1.split("_")[4][2:4]),
        second=int(basef1.split("_")[4][4:6]),
    )

    secondary["begintime"] = datetime(
        year=int(basef2.split("_")[3][6:10]),
        month=int(basef2.split("_")[3][10:12]),
        day=int(basef2.split("_")[3][12:14]),
        hour=int(basef2.split("_")[4][0:2]),
        minute=int(basef2.split("_")[4][2:4]),
        second=int(basef2.split("_")[4][4:6]),
    )

    print(
        "Reference:\n  Name:          "
        + reference["filename"]
        + "\n  Begin time:    "
        + reference["begintime"].isoformat()
        + "\n  Platform:      "
        + reference["platform"]
        + "\n  File location: "
        + reference["fileloc"]
    )
    print(
        "\nSecondary:\n  Name:          "
        + secondary["filename"]
        + "\n  Begin time:    "
        + secondary["begintime"].isoformat()
        + "\n  Platform:      "
        + secondary["platform"]
        + "\n  File location: "
        + secondary["fileloc"]
    )

    timedelta = secondary["begintime"] - reference["begintime"]

    print(
        "\nPair information:\n  Baseline:      "
        + str(timedelta)
        + "\n  Track:         "
        + str(reference["track"])
        + "\n  Orbits:        "
        + str(reference["orbit"])
        + ", "
        + str(secondary["orbit"])
    )

    reference["begintime"] = reference["begintime"].isoformat()
    secondary["begintime"] = secondary["begintime"].isoformat()

    logger.addFrameMetadata("reference", reference)
    logger.addFrameMetadata("secondary", secondary)

    # reference
    if reference["platform"] == "ERS-1":
        reference_orbitloc = "/home/data/orbits/ODR/ERS1"

        reference_orbit = f"""<property name="ORBIT_TYPE">
            <value>ODR</value>
        </property>
        <property name="ORBIT_DIRECTORY">
            <value>{reference_orbitloc}</value>
        </property>"""

        reference_file = f"""<property name="IMAGEFILE">
            <value>{reference['fileloc']}</value>
        </property>"""

        sensor_name = "ERS_EnviSAT"
    elif reference["platform"] == "ERS-2":
        reference_orbitloc = "/home/data/orbits/ODR/ERS2"

        reference_orbit = f"""<property name="ORBIT_TYPE">
            <value>ODR</value>
        </property>
        <property name="ORBIT_DIRECTORY">
            <value>{reference_orbitloc}</value>
        </propertcdy>"""

        reference_file = f"""<property name="IMAGEFILE">
            <value>{reference['fileloc']}</value>
        </property>"""

        sensor_name = "ERS_EnviSAT"
    elif reference["platform"] == "Envisat":
        reference_file = f"""<property name="IMAGEFILE">
            <value>{reference['fileloc']}</value>
        </property>"""

        reference_orbit = f"""<property name="INSTRUMENTFILE">
            <value>"/home/yadevries/sarpixeltracking/data/instrumentfiles/ASA_INS_AXVIEC20061220_105425_20030211_000000_20071231_000000"</value> 
        </property>
        <property name="ORBITFILE">
            <value>"/home/yadevries/sarpixeltracking/data/orbitfiles/DOR_VOR_AXVF-P20201101_200200_20031115_215528_20031117_002328"</value> 
        </property>"""

        sensor_name = "EnviSAT"

    print("\nreference.xml:")

    reference_xml = f"""
    <component name="Reference">
        {reference_file}
        <property name="OUTPUT">reference</property>
        {reference_orbit}
    </component>"""

    print(reference_xml)

    f = open("reference.xml", "w")
    f.write(reference_xml)
    f.close()

    # secondary
    if secondary["platform"] == "ERS-1":
        secondary_orbitloc = "/home/data/orbits/ODR/ERS1"

        secondary_fileloc = f"""<property name="IMAGEFILE">
            <value>{secondary['fileloc']}</value>
        </property>"""
    elif secondary["platform"] == "ERS-2":
        secondary_orbitloc = "/home/data/orbits/ODR/ERS2"

        secondary_fileloc = f"""<property name="IMAGEFILE">
            <value>{secondary['fileloc']}</value>
        </property>"""
    elif secondary["platform"] == "Envisat":
        secondary_file = f"""<property name="IMAGEFILE">
            <value>{secondary['fileloc']}</value>
        </property>"""

        secondary_orbit = f"""<property name="INSTRUMENTFILE">
            <value>"/home/yadevries/sarpixeltracking/data/instrumentfiles/ASA_INS_AXVIEC20061220_105425_20030211_000000_20071231_000000"</value> 
        </property>
        <property name="ORBITFILE">
            <value>"/home/yadevries/sarpixeltracking/data/orbitfiles/DOR_VOR_AXVF-P20201101_203000_20040124_215528_20040126_002328"</value> 
        </property>"""

    print("\nsecondary.xml:")

    secondary_xml = f"""
    <component name="Secondary">
        {secondary_file}
        <property name="OUTPUT">secondary</property>
        {secondary_orbit}
    </component>"""

    print(secondary_xml)

    f = open("secondary.xml", "w")
    f.write(secondary_xml)
    f.close()

    # stripmapApp.xml
    print("\nstripmapApp.xml:")

    # DEM_loc = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem"  # "/home/yad2/DEM/IslandsDEMv1.0_2x2m_zmasl_isn93_SouthMerge.tif"

    stripmapApp_xml = f"""
    <stripmapApp>
        <component name="insar">
            <property  name="Sensor name">{sensor_name}</property>
            <component name="reference">
                <catalog>reference.xml</catalog>
            </component>
            <component name="secondary">
                <catalog>secondary.xml</catalog>
            </component>
            <property name="demFilename">{demfile}</property>
            <property name="do denseoffsets">True</property>
            <property name="geocode list">['filt_topophase.flat', 'los.rdr', 'topophase.cor', 'phsig.cor', 'reference_slc/reference.slc', 'secondary_slc/secondary.slc']</property>
        </component>
    </stripmapApp>"""

    print(stripmapApp_xml)

    f = open("stripmapApp.xml", "w")
    f.write(stripmapApp_xml)
    f.close()

    # stripmapApp.xml
    print("\ndense.xml:")

    # DEM_loc = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem"  # "/home/yad2/DEM/IslandsDEMv1.0_2x2m_zmasl_isn93_SouthMerge.tif"

    dense_xml = f"""
    <stripmapAppDenseAmpcor>
        <component name="dense">
            <property name="Ampcor window width">64</property>
            <property name="Ampcor window height">256</property>
            <property name="Ampcor search window width">10</property>
            <property name="Ampcor search window height">40</property>
            <property name="Ampcor skip width">128</property>
            <property name="Ampcor skip height">32</property>
        </component>
    </stripmapAppDenseAmpcor>"""

    print(dense_xml)

    f = open("dense.xml", "w")
    f.write(dense_xml)
    f.close()


def generate_files_TSX_format(logger, file1, file2, demfile):
    
    import os
    from datetime import datetime

    reference = {}
    secondary = {}

    reference["fileloc"] = file1
    secondary["fileloc"] = file2

    reference["filename"] = os.path.basename(file1)
    secondary["filename"] = os.path.basename(file2)

    basef1, extf1 = os.path.splitext(os.path.basename(file1))

    basef2, extf2 = os.path.splitext(os.path.basename(file2))
    
    reference["platform"] = "TSX"
    secondary["platform"] = "TSX"
        
    reference["begintime"] = datetime(
        year=int(basef1.split("_")[-2][:4]),
        month=int(basef1.split("_")[-2][4:6]),
        day=int(basef1.split("_")[-2][6:8]),
        hour=int(basef1.split("_")[-2][9:11]),
        minute=int(basef1.split("_")[-2][11:13]),
        second=int(basef1.split("_")[-2][13:15]),
    )
    
    secondary["begintime"] = datetime(
        year=int(basef2.split("_")[-2][:4]),
        month=int(basef2.split("_")[-2][4:6]),
        day=int(basef2.split("_")[-2][6:8]),
        hour=int(basef2.split("_")[-2][9:11]),
        minute=int(basef2.split("_")[-2][11:13]),
        second=int(basef2.split("_")[-2][13:15]),
    )
    
    print(
        "Reference:\n  Name:          "
        + reference["filename"]
        + "\n  Begin time:    "
        + reference["begintime"].isoformat()
        + "\n  Platform:      "
        + reference["platform"]
        + "\n  File location: "
        + reference["fileloc"]
    )
    print(
        "\nSecondary:\n  Name:          "
        + secondary["filename"]
        + "\n  Begin time:    "
        + secondary["begintime"].isoformat()
        + "\n  Platform:      "
        + secondary["platform"]
        + "\n  File location: "
        + secondary["fileloc"]
    )

    timedelta = secondary["begintime"] - reference["begintime"]

    print(
        "\nPair information:\n  Baseline:      "
        + str(timedelta)
    )
    
    reference["begintime"] = reference["begintime"].isoformat()
    secondary["begintime"] = secondary["begintime"].isoformat()

    logger.addFrameMetadata("reference", reference)
    logger.addFrameMetadata("secondary", secondary)
    
    
    print("\nreference.xml:")

    reference_xml = f"""
    <component name="Reference">
        <property name="XML">
            {reference["fileloc"]}
        </property>
        <property name="OUTPUT">reference</property>
    </component>"""

    print(reference_xml)

    f = open("reference.xml", "w")
    f.write(reference_xml)
    f.close()

    print("\nsecondary.xml:")

    secondary_xml = f"""
    <component name="Secondary">
        <property name="XML">
            {secondary["fileloc"]}
        </property>
        <property name="OUTPUT">secondary</property>
    </component>"""

    print(secondary_xml)

    f = open("secondary.xml", "w")
    f.write(secondary_xml)
    f.close()

    # stripmapApp.xml
    print("\nstripmapApp.xml:")

    # DEM_loc = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem"  # "/home/yad2/DEM/IslandsDEMv1.0_2x2m_zmasl_isn93_SouthMerge.tif"

    stripmapApp_xml = f"""
    <stripmapApp>
        <component name="insar">
            <property  name="Sensor name">TerraSARX</property>
            <component name="reference">
                <catalog>reference.xml</catalog>
            </component>
            <component name="secondary">
                <catalog>secondary.xml</catalog>
            </component>
            <property name="demFilename">{demfile}</property>
            <property name="do denseoffsets">True</property>
            <property name="regionOfInterest">[63.699855,63.583704,-19.476357,-19.205132]</property>
            <property name="geocode list">['interferogram/filt_topophase.flat', 'interferogram/los.rdr', 'interferogram/topophase.cor', 'interferogram/phsig.cor', 'reference_slc/reference.slc', 'coregisteredSlc/refined_coreg.slc']</property>
        </component>
    </stripmapApp>"""

    print(stripmapApp_xml)

    f = open("stripmapApp.xml", "w")
    f.write(stripmapApp_xml)
    f.close()

    # # stripmapApp.xml
    # print("\ndense.xml:")

    # # DEM_loc = "/home/data/DEM/LMI/ArcticDEM/v1/Iceland_10m.dem"  # "/home/yad2/DEM/IslandsDEMv1.0_2x2m_zmasl_isn93_SouthMerge.tif"

    # dense_xml = f"""
    # <stripmapAppDenseAmpcor>
    #     <component name="dense">
    #         <property name="Ampcor window width">64</property>
    #         <property name="Ampcor window height">256</property>
    #         <!--<property name="Ampcor search window width">10</property>-->
    #         <!--<property name="Ampcor search window height">40</property>-->
    #         <property name="Ampcor skip width">128</property>
    #         <property name="Ampcor skip height">32</property>
    #     </component>
    # </stripmapAppDenseAmpcor>"""

    # print(dense_xml)

    # f = open("dense.xml", "w")
    # f.write(dense_xml)
    # f.close()