"""
This is just a simple test script to verify that 
I can properly read mapping files:
"""
import sys, os, glob
from cpl import map

# Get some files to read:
files = glob.glob("data/*.nc")

maps = []
for file in files:
    mymap = map.Map()
    mymap.read(file)
    maps.append(mymap)




