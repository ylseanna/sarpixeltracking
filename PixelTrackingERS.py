#!/usr/bin/env python

'''
Base function for pixel tracking for this project

Requires a file called MetaData.db containting the file locations of all files used in this project
'''

### imports

import os, glob
import sqlite3 #database
from datetime import datetime #timemanagement


### utils

def timestamp(string):
    return datetime.strptime(string, "%Y-%m-%dT%H:%M:%S")


### INPUT

print("\n*** ESR 1&2 - Pixel Tracking Pipeline ***\n")

program_start = datetime.now()

print(program_start)

print("\n - Initialising...\n")

'''
EVENTUALLY MAKE IT A COMMAND LINE ARGUMENT, for now hardcoded
'''


file1 = "ER01_SAR_IM__0P_19920524T123500_19920524T123517_UPA_04477_0000.CEOS.tar.gz"
file2 = "ER01_SAR_IM__0P_19930718T123508_19930718T123525_UPA_10489_0000.CEOS.tar.gz"

print("Files selected:\n"+file1+"\n"+file2)


### QUERY FILES and ASSIGN REF AND SEC

print("\n - Querying files...")

def query(query_str):
    conn = sqlite3.connect("MetaData.db")

    cursor = conn.execute(query_str)

    return {'results':
            [dict(zip([column[0] for column in cursor.description], row))
             for row in cursor.fetchall()]}

    conn.close()

files = query("SELECT * from frames WHERE filename = '{}' OR filename = '{}'".format(file1, file2))

file1, file2 = files['results']

print("\nFiles succesfully queried:\n")

timedelta = timestamp(file1['begintime']) - timestamp(file2['begintime'])

if timedelta.total_seconds() < 0:
    reference = file1
    secondary = file2
else:
    reference = file2
    secondary = file1


print("Reference:\n  Name:          "+reference['filename']+"\n  Begin time:    "+reference['begintime']+"\n  End time:      "+reference['endtime']+"\n  Sensor:        "+reference['sensor']+"\n  File location: "+reference['fileloc'])
print("\nSecondary:\n  Name:          "+secondary['filename']+"\n  Begin time:    "+secondary['begintime']+"\n  End time:      "+secondary['endtime']+"\n  Sensor:        "+secondary['sensor']+"\n  File location: "+secondary['fileloc'])

timedelta = timestamp(secondary["begintime"]) - timestamp(reference["begintime"])

print("\nPair information:\n  Baseline:      "+str(timedelta)+"\n  Track:         "+str(reference['track'])+"\n  Frame:         "+str(reference['frame'])+"")

### CREATE XML files and processing folders

### /home/data/orbits/ODR/ERS1
### /home/data/orbits/ODR/ERS2

# Folders:

print("\n - Generating folder structure...")

#make sure ref_dir exists and is empty
if os.path.exists("./reference") == False:
    os.mkdir("reference")
else: 
    files = glob.glob("./reference/*")
    for f in files:
        os.remove(f)

#make sure sec_dir exists and is empty
if os.path.exists("./secondary") == False:
    os.mkdir("secondary")
else: 
    files = glob.glob("./secondary/*")
    for f in files:
        os.remove(f)

# Unpack images:

print("\n - Unpacking images...")

print("\nReference:\n")

os.system(f"tar -zvxf {reference['fileloc']} --directory ./reference")

print("\nSecondary:\n")

os.system(f"tar -zvxf {secondary['fileloc']} --directory ./secondary")

# Generate XML-files

print("\n - Generating XML-files...")


#reference
if reference['sensor'] == 'ERS 1':
    reference_orbitloc = "/home/data/orbits/ODR/ERS1"
elif reference["sensor"] == 'ERS 2':
    reference_orbitloc = "/home/data/orbits/ODR/ERS2"

print("\nreference.xml:")


reference_xml = f'''<component name="Reference">
    <property name="IMAGEFILE">
        ./reference/DAT_01.001
    </property>
    <property name="LEADERFILE">
        ./reference/LEA_01.001
    </property>
    <property name="OUTPUT">reference.raw</property>
    <property name="ORBIT_TYPE">
        <value>ODR</value>
    </property>
    <property name="ORBIT_DIRECTORY">
        <value>{reference_orbitloc}</value>
    </property>
</component>'''

print(reference_xml)

f = open("reference.xml", "w")
f.write(reference_xml)
f.close()

#secondary
if secondary['sensor'] == 'ERS 1':
    secondary_orbitloc = "/home/data/orbits/ODR/ERS1"
elif secondary["sensor"] == 'ERS 2':
    secondary_orbitloc = "/home/data/orbits/ODR/ERS2"

print("\nsecondary.xml:")

secondary_xml = f'''<component name="Secondary">
    <property name="IMAGEFILE">
        ./secondary/DAT_01.001
    </property>
    <property name="LEADERFILE">
        ./secondary/LEA_01.001
    </property>
    <property name="OUTPUT">secondary.raw</property>
    <property name="ORBIT_TYPE">
        <value>ODR</value>
    </property>
    <property name="ORBIT_DIRECTORY">
        <value>{secondary_orbitloc}</value>
    </property>
</component>'''

print(secondary_xml)

f = open("secondary.xml", "w")
f.write(secondary_xml)
f.close()



#stripmapApp.xml
print("\nstripmapApp.xml:")

stripmapApp_xml = f'''<stripmapApp>
    <component name="insar">
        <property  name="Sensor name">ERS</property>
        <component name="reference">
            <catalog>reference.xml</catalog>
        </component>
        <component name="secondary">
            <catalog>reference.xml</catalog>
        </component>
    </component>
</stripmapApp>'''

print(stripmapApp_xml)

f = open("stripmapApp.xml", "w")
f.write(stripmapApp_xml)
f.close()





# Finish

program_end = datetime.now()

print(f"\nProcessing finished at {program_end}")
print(f"Completed in {program_end - program_start}")


