#!/usr/bin/python
# Ryan Valenza
# 2014-10-28
# Mod 2016-08-08 MWistey: fix regexp, replace eval->enval
# This script is the start of a project to create a series of
# modules to interpret VASP outfiles.  The end goal is to be able
# to quickly create plots and format data for external use.

import sys
import re
import math

try:
	eigenval = open("EIGENVAL","r")
except IOError:
	sys.exit("Could not open EIGENVAL")

eigenval.readline() # N/A
eigenval.readline() # N/A
eigenval.readline() # N/A
eigenval.readline() # Cartesian/Direct

# System name w/ stripped newline character
name = eigenval.readline().rstrip() 
print "# System:        " + name

first = eigenval.readline() # Possibly interesting information
(nelect,nkpts,nbands) = first.split()
print "# of electrons:  " + nelect
print "# of k-points:   " + nkpts
print "# of bands:      " + nbands

# Regular expressions are used to distinguish between a k-point and an eigenvalue
regexs = {
'kpt': "\s+(.\d+\.\d+E[+-]\d+)\s+(.\d+\.\d+E[+-]\d+)\s+(.\d+\.\d+E[+-]\d+)\s+.*",
'enval': "\s+(\d+)\s+(-?.\d+\.\d+)" }

kpts = []
bands = []
for i in range(int(nbands)):
	bands.append([]) # lists are immutable
j = 0 # mark band number

print "Finding k-point, energy pairs for each band..."
for line in eigenval:
	kpt = re.match(regexs['kpt'], line)
	enval = re.match(regexs['enval'], line)

	if kpt != None:
		(kx,ky,kz) = kpt.groups()
		k = math.sqrt(float(kx)**2 + float(ky)**2 + float(kz)**2)
		k = k
		
	if enval != None:
		e = float(enval.groups(0)[1])
		bands[j%int(nbands)].append([k,e])
		j += 1

# As of right now, the data being is being organized into a format easily read
# into Mathematica.  The format can be changed to fit future problems.
# Ideas - use matplotlib to generate a bandstructure plot
#       - numpy arrays
print "Creating file bands.txt..."
out = open("bands.txt", "w")
for i in range(int(nbands)):
	for j in range(int(nkpts)):
		out.write("%d %.4f %.4f\n"%(i,bands[i][j][0],bands[i][j][1]))


