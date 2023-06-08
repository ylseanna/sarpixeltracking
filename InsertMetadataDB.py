#!/usr/bin/env python

'''
Just a simple script to import useable image file locations for the Tungnarkvislarjokull landslide project into the metadata database, currently set up to import files from the GEANT project in the folder ERS_new on the data disk.
'''

import os
import sqlite3

#variables
root = "/home/data/ERS_new/SLC/T324"
track = 341
frames = [
	2277,
	2295,
	2313
]


def TransferTime(string):
	year 	= string[:4]
	month 	= string[4:6]
	day 	= string[6:8]

	hour 	= string[9:11]
	minute 	= string[11:13]
	second	= string[13:15]

	return year+"-"+month+"-"+day+"T"+hour+":"+minute+":"+second


for path, subdirs, files in os.walk(root):
	for i, f in enumerate(files):
		if f[:4] == "ER01":
			sensor = "ERS 1"
		elif f[:4] == "ER02":
			sensor = "ERS 2"

		datalevel 	= f[13:15]
		begintime 	= TransferTime(f[16:31])
		endtime 	= TransferTime(f[32:47])
		orbit 		= int(f[52:57])
                
		filepath	= os.path.join(path, f)

		# print(i, f, begintime, endtime, track, frames[i], sensor, orbit, datalevel, filepath)


		#database insertion
		conn = sqlite3.connect("MetaData.db")


		keysvalues 	= {
            "filename"  :  f,  
            "begintime" :  begintime,
            "endtime"   :  endtime,
            "track"     :  track,
            "frame"     :  frames[i],
            "sensor"    :  sensor,     
            "orbit"     :  orbit,
            "datalevel" :  datalevel,
            "fileloc"   :  filepath
        }

		#prepare statement
		columns 	= "("
		values  	= "("

		for key, value in keysvalues.items():
			columns += key+","

			if type(value) == float:
				values  += str(value)+","
			else:
				values  += "'"+str(value)+"',"


		columns = columns[:-1] + ")"
		values  = values[:-1] + ")"
	
		statement = "INSERT INTO frames {} VALUES {};".format(columns, values)

		print(statement)

		conn.execute(statement)

		conn.commit()   

		conn.close()


