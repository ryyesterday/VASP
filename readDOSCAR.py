#!/usr/bin/python
# Ryan Valenza
# 2014-11-06
# This program is the start of a project to create a series of
# modules to interpret VASP outfiles.  The end goal is to be able
# to quickly create plots and format data for external use.

# Note:  It has been assumed that the VASP calculation was non
#		 collinear (non spin-polarized).  This may be implemented
#		 by setting ISPIN = 1 in INCAR.  Therefore, the default 
#		 format is [energy] [DOS-total] [DOS-integrated]

import matplotlib.pyplot as plt

try:
	doscar = open("DOSCAR","r")
except IOERROR:
	sys.exit("Could not open DOSCAR")

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

plt.figure(1)			# The first figure
plt.subplot(2,1,1)		# The first subplot in the first figure
plt.plot(es,tds,'b-')
plt.ylabel('Total DOS')

plt.subplot(2,1,2)		# The second subplot in the first figure
plt.plot(es,ids,'r-')
plt.ylabel('Integrated DOS')

plt.xlabel(r'$E - E_f$ (eV)')
plt.subplot(2,1,1) 		# Make first subplot current (for title)
plt.title(name + " Density of States")
plt.show()