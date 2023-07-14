#!/usr/bin/env python

'''
Just a simple script to import useable TSX image file locations for the Tungnarkvislarjokull landslide project into the metadata database, currently set up to import files from the GEANT project in the folder ERS_new on the data disk.
'''

import os
import sqlite3
import glob

#variables
root = "/home/data/TSX/South/T41_s011_6345_Katla"
track = 41
frame = 6345


def TransferTime(string):
	year 	= string[:4]
	month 	= string[4:6]
	day 	= string[6:8]

	hour 	= string[9:11]
	minute 	= string[11:13]
	second	= string[13:15]

	return year+"-"+month+"-"+day+"T"+hour+":"+minute+":"+second

folders = [f.path for f in os.scandir(root) if f.is_dir()]


for folder in folders:
	result1 = glob.glob(os.path.join(folder, 'TSX-1.SAR.L1B', '*'))

	subfolder = result1[0].split('/')[-1]

	filepath = glob.glob(os.path.join(folder, 'TSX-1.SAR.L1B', subfolder, '*.xml'))[0]

	f = filepath.split('/')[-1]

	print(f)

	sensor = 'TSX1'

	timestrings = f[28:-4]

	datalevel 	= 'na'
	begintime 	= TransferTime(timestrings.split('_')[0])
	endtime 	= TransferTime(timestrings.split('_')[1])
	orbit 		= 'na'
			

	# print(i, f, begintime, endtime, track, frames[i], sensor, orbit, datalevel, filepath)


	#database insertion
	conn = sqlite3.connect("MetaData.db")


	keysvalues 	= {
		"filename"  :  f,  
		"begintime" :  begintime,
		"endtime"   :  endtime,
		"track"     :  track,
		"frame"     :  frame,
		"sensor"    :  sensor,     
		"orbit"     :  orbit,
		"datalevel" :  datalevel,
		"fileloc"   :  filepath
	}

	# print(keysvalues)

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


