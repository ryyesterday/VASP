#!/usr/bin/python
# Ryan Valenza
# 2014-11-06
# This script is the start of a project to create a series of
# modules to interpret VASP outfiles.  The end goal is to be able
# to quickly create plots and format data for external use.

# Note:  It has been assumed that the VASP calculation was non
#		 collinear (non spin-polarized).  This may be implemented
#		 by setting ISPIN = 1 in INCAR.  Therefore, the default 
#		 format is [energy] [DOS-total] [DOS-integrated]

import matplotlib.pyplot as plt
import sys

def readDOSCAR(filename):

	doscar = open(filename,"r")
	filename = filename.split("\\")[-1]
	# e.g. /home/ryval/DOSCAR, grab DOSCAR
	print "# Filename:  " + filename
    
	# Skip file header
	doscar.readline() # N/A
	doscar.readline() # N/A
	doscar.readline() # N/A
	doscar.readline() # Cartesian/Direct
	
	# System name w/ stripped newline character
	name = doscar.readline().rstrip() 
	print "# System:	" + name 

	first = doscar.readline()
	(emax, emin, ndos, efermi, weight) = first.split()
	print "# Max Energy:    %.4f" % float(emax) + " eV"
	print "# Min Energy:    %.4f" % float(emin) + " eV"
	print "# Fermi Energy:  %.4f" % float(efermi) + " eV"
	print "# " + ndos + " points calculated"
	print # Buffer

	es = []
	tds = []
	ids = []

	# Loop through the rest of the file and grab DOS - energy pairs
	for i in range(int(ndos)):
		line = doscar.readline()
		(energy, tdos, idos) = line.split()
		es.append(float(energy)-float(efermi))
		tds.append(float(tdos))
		ids.append(float(idos))

	plt.subplot(2,1,1)		# The first subplot in the first figure
	plt.plot(es,tds,'-')
	plt.ylabel('Total DOS')

	plt.subplot(2,1,2)		# The second subplot in the first figure
	plt.plot(es,ids,'-')
	plt.ylabel('Integrated DOS')

	plt.xlabel(r'$E - E_f$ (eV)')
	plt.subplot(2,1,1) 		# Make first subplot current (for title)
	plt.title(name + " Density of States")
	
	
def main():

	"""
readDOSCAR.py is used to generate plots of the total 
and integrated density of states calculated using VASP.  
it can be passed a list of filenames via CLI or a file 
that contains a list of filenames, using the -f option.

readDOSCAR.py DOSCAR1 DOSCAR2 DOSCAR3 ...				
readDOSCAR.py -f [DOS-file]				 

Example DOS-file:						
/path/to/DOSCAR1
/path/to/DOSCAR2
/path/to/DOSCAR3 
	"""

	print # Buffer
	plt.figure('Density of States', figsize = (10,10)) # Create figure 1
	
	if sys.argv[1] == '-f':
		filelist = open(sys.argv[2],"r")
		for dosfile in filelist:
			readDOSCAR(dosfile.rstrip())
	else:
		for cl_arg in sys.argv:
			if cl_arg != sys.argv[0]: # Skip readDOSCAR.py from CL
				readDOSCAR(cl_arg)
			
	plt.show() # Display the figure

	
# Run the script, if it breaks, print a usage statement
try:
	main()
except IOError: # File error
	print "ERROR - Unable to open DOSCAR file"
	print main.__doc__
except:	# Anything else
	print "ERROR - Are you using readDOSCAR.py correctly?"
	print main.__doc__
	