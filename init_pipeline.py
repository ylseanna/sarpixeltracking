
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
    os.system("rm *.xml")


    print("Files selected:\n" + file1 + "\n" + file2)

    ### QUERY FILES and ASSIGN REF AND SEC

    print("\n - Querying files...")
    
    # if file1.endswith(".CEOS.tar.gz"):
        # generate_files_ERS_CEOS(logger, file1, file2, demfile) # NEEDS TO BE ADJUSTED TO NOT USE THE DATABASE
    if file1.endswith((".E1", ".E2", ".N1")):
        generate_files_Envisat_format(logger, file1, file2, demfile)



    
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
        
    reference["track"], secondary["track"] = (int(basef1.split("_")[6]), int(basef2.split("_")[6]))
    
    reference["orbit"], secondary["orbit"] = (int(basef1.split("_")[7]), int(basef2.split("_")[7]))
    
    reference["begintime"] = datetime(
        year = int(basef1.split("_")[3][6:10]),
        month = int(basef1.split("_")[3][10:12]),
        day = int(basef1.split("_")[3][12:14]),
        hour = int(basef1.split("_")[4][0:2]),
        minute = int(basef1.split("_")[4][2:4]),
        second = int(basef1.split("_")[4][4:6])
    )
    
    secondary["begintime"] = datetime(
        year = int(basef2.split("_")[3][6:10]),
        month = int(basef2.split("_")[3][10:12]),
        day = int(basef2.split("_")[3][12:14]),
        hour = int(basef2.split("_")[4][0:2]),
        minute = int(basef2.split("_")[4][2:4]),
        second = int(basef2.split("_")[4][4:6])
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
    elif reference["platform"] == "ERS-2":
        reference_orbitloc = "/home/data/orbits/ODR/ERS2"

    print("\nreference.xml:")

    reference_xml = f"""
    <component name="Reference">
        <property name="IMAGEFILE">
            <value>{reference['fileloc']}</value>
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
            <value>{secondary['fileloc']}</value>
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
            <property  name="Sensor name">ERS_EnviSAT</property>
            <component name="reference">
                <catalog>reference.xml</catalog>
            </component>
            <component name="secondary">
                <catalog>secondary.xml</catalog>
            </component>
            <property name="demFilename">{demfile}</property>
            <property name="do denseoffsets">True</property>
            <!--<property name="regionOfInterest">[63.615914,63.697878,-19.500389,-19.240837]</property>-->
        </component>
    </stripmapApp>"""

    print(stripmapApp_xml)

    f = open("stripmapApp.xml", "w")
    f.write(stripmapApp_xml)
    f.close()
  
    
    
