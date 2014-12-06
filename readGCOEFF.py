#!/usr/bin/python
# Ryan Valenza
# 2014-12-05
# This script is the start of a project to create a series of
# modules to interpret VASP outfiles.  The end goal is to be able
# to quickly create plots and format data for external use.

# Imports
import matplotlib.pyplot as plt
import sys
import time

# Global Variables
OCCTOL = 10 ** (-6) # occupation tolerance
OUTFILE = open("C:\Users\Ryan\Desktop\cleaned_GCOEFF.txt",'w')

# The readGCOEFF class
# Must initialize with a GCOEFF filename
class readGCOEFF:

	# Initialize the reader object
	def __init__(self, filename):
	
		# open the file, 1 stands for buffer
		self.gcoeff = open(filename,'r', 1)
		# for looping & such
		self.ispin = 0
		self.nkpts = 0
		self.nbands = 0
		self.A = []
		self.B = []
		self.data = []
	
	# Read the header
	def readHeader(self):

		self.ispin = int(self.gcoeff.readline())
		self.nkpts = int(self.gcoeff.readline())
		self.nbands = int(self.gcoeff.readline())
		# Lattice vector
		self.A.append([float(i) for i in self.gcoeff.readline().split()])
		self.A.append([float(i) for i in self.gcoeff.readline().split()])
		self.A.append([float(i) for i in self.gcoeff.readline().split()])
		# Reciprocal vector
		self.B.append([float(i) for i in self.gcoeff.readline().split()])
		self.B.append([float(i) for i in self.gcoeff.readline().split()])
		self.B.append([float(i) for i in self.gcoeff.readline().split()])
	
	# W/in an unoccupied band, skip all PWs
	def __skipPWs(self):
		
		nps = int(self.npws)
		while nps >= 1:
			self.gcoeff.readline()
			nps -= 1
	
	# W/in a band, read all PWs
	def __readPWs(self):
	
		nps = int(self.npws)
		while nps >= 1:
		
			gline = self.gcoeff.readline()
			(g1, g2, g3, d1, real, d2, im, d3) = gline.split()
			
			# save data to object for future
			self.data.append([self.n,			# band number
				self.k1,self.k2,self.k3,	    # k-pt
				g1,g2,g3,						# G-pt
				real,im,						# PW coefficient
				self.energy])					# energy
				
			# write data to outfile for mathematica import
			OUTFILE.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n" % (self.n,
				self.k1,self.k2,self.k3,
				g1,g2,g3,
				real,im,
				self.energy))
			
			nps -= 1
			
	# W/in a kpt, read all band data
	def __readBand(self):
	
		nbds = self.nbands
		while nbds >= 1:
		
			nline = self.gcoeff.readline()
			(self.n, self.npws) = nline.split()
			
			# check occupancy
			occline = self.gcoeff.readline()
			self.energy = occline.split()[1]
			self.occ = float(occline.split()[5])
			if self.occ >= OCCTOL:
				self.__readPWs()
			else:
				self.__skipPWs()
				
			nbds -= 1
	
	# Take k1,k2,k3 and multiple by B
	def __getK(self):
	
		k1 = [float(self.k1) * i for i in self.B[0]]
		k2 = [float(self.k2) * i for i in self.B[1]]
		k3 = [float(self.k3) * i for i in self.B[2]]
		k = [x+y+z for x,y,z in zip(k1,k2,k3)]
		
		# reform k
		self.k1 = k[0]
		self.k2 = k[1]
		self.k3 = k[2]
		
	# Read the entire file
	def readFile(self):

		self.readHeader()

		nks = self.nkpts
		while nks >= 1:
		
			k = self.gcoeff.readline()
			(self.k1, self.k2, self.k3) = k.split()
			self.__getK() # form the actual k-vector from B
			self.__readBand()
			nks -= 1
			

# use the class			

def main():

	start = time.time()		

	filename = sys.argv[1]
	obj = readGCOEFF(filename)
	obj.readFile()

	end = time.time()

	print "Time elapsed:  %4f"%(end-start)
	
	
main()

