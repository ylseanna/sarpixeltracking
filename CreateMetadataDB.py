#!/usr/bin/env python

'''
This function generates the metadata (frames and most importantly frame location) database for the pixel tracking pipeline.

This function can be fully run to completely reset the database.
'''

#create database
import sqlite3
dbName = 'MetaData.db'

try:
  conn = sqlite3.connect(dbName)
  cursor = conn.cursor()
  print("Database created!")

except Exception as e:
  print("Something bad happened: ", e)
  if conn:
    conn.close()

# Create table for frames
create_query = '''
    CREATE TABLE IF NOT EXISTS frames(
    id         INTEGER PRIMARY KEY,
    filename   TEXT NOT NULL,
    begintime  TEXT NOT NULL,
    endtime    TEXT NOT NULL,
    track      FLOAT NOT NULL,
    frame      FLOAT NOT NULL,
    sensor     TEXT NOT NULL,
    orbit      FLOAT NOT NULL,
    datalevel  TEXT NOT NULL,
    fileloc    TEXT NOT NULL);
'''

cursor.execute(create_query)
print("Table for frames created!")


# Create table for image pairs
create_query = '''
    CREATE TABLE IF NOT EXISTS pairs(
    id             INTEGER PRIMARY KEY,
    referencename  TEXT NOT NULL,
    referenceid    INTEGER NOT NULL,
    secondaryname  TEXT NOT NULL,
    secondaryid    INTEGER NOT NULL,
    referencetime  TEXT NOT NULL,
    secondarytime  TEXT NOT NULL,
    track          FLOAT NOT NULL,
    frame          FLOAT NOT NULL,
    sensor         TEXT NOT NULL,
    fileloc        TEXT NOT NULL);
'''

cursor.execute(create_query)
print("Table for image pairs created!")

cursor.close()