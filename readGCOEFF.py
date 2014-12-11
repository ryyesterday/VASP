#!/usr/bin/python
# Ryan Valenza
# 2014-12-05
# This script is the start of a project to create a series of
# modules to interpret VASP outfiles.  The end goal is to be able
# to quickly create plots and format data for external use.

# Imports
from numpy import *
import sys
from sets import Set
import time

# Global Variables
OCCTOL = 10 ** (-6) # occupation tolerance
OUTFILE = "C:\Users\Ryan\Desktop\cleaned_GCOEFF.txt"
SFFILE = "C:\Users\Ryan\Desktop\SFactors.txt"

# The readGCOEFF class
# Must initialize with a GCOEFF filename
class readGCOEFF:


	# Initialize the reader object
	def __init__(self, filename):
	
		# open the file, 1 stands for buffer
		self.gcoeff = open(filename,'r', 1)
		
		# for looping & such
		self.start = True
		self.ispin = 0
		self.nkpts = 0
		self.nbands = 0
		
		# data storage
		self.A = []
		self.B = []
		self.data = []
		self.big_data = []
		
		# dictionary for h,k,l lookup, constant time
		self.coefs = {}
		self.big_coefs = []
		self.G = Set()
		
		# prev, used to remove querying for n & k
		self.prevN = 1
		self.prevKx = 0
		self.prevKy = 0
		self.prevKz = 0
		
		# outfile and checks
		self.write_out = True
		self.cartK = False
		self.read = False
		self.readH = False
	
	
	# given a segmented data set, find a coefficient given G
	# segmented = unique n & k; input must be floats
	def __findCoeff(self, seg, gx, gy, gz):
		
		gc = seg[logical_and(seg[:,4] == gx,
			logical_and(seg[:,5] == gy,
			seg[:,6] == gz))]
			
		# check if NULL result
		if gc.shape[0] == 0:
			return [0, 0]
		else:
			# return [real, imag] PW coefficient
			return [gc[0][7],gc[0][8]]
		
		
	# get structure factor, S(Q)
	def __getSF(self, qx, qy, qz):
	
		qx = int(qx)
		qy = int(qy)
		qz = int(qz)
		
		sf = [0,0]
		
		for coef in self.big_coefs:
			for g in self.G:
			
				gx = g[0]
				gy = g[1]
				gz = g[2]
				
				try:
					c1 = coef['[%i,%i,%i]'%(gx,gy,gz)]
					c2 = coef['[%i,%i,%i]'%(gx+qx,gy+qy,gz+qz)]
				except KeyError:
					continue
				
				# form real & imaginary part
				re = c1[0]*c2[0] + c1[1]*c2[1]
				im = c1[1]*c2[0] - c1[0]*c2[1]
				sf[0] += re
				sf[1] += im
				
		return sf
		
	# Read the header
	def readHeader(self):

		print "Reading header..."
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
		
		nps = self.npws
		while nps >= 1:
			self.gcoeff.readline()
			nps -= 1
	
	
	# W/in a band, read all PWs
	def __readPWs(self):
	
		nps = self.npws
		while nps >= 1:
		
			gline = self.gcoeff.readline()
			(g1, g2, g3, d1, re, d2, im, d3) = gline.split()
			re = float(re)
			im = float(im)
			g1 = int(g1)
			g2 = int(g2)
			g3 = int(g3)
			
			self.G.add((g1,g2,g3))  # add to set of all G-values
			
			# if current = prev, append to existing dataset
			if (self.prevKx == self.k1 and 
				self.prevKy == self.k2 and
				self.prevKz == self.k3 and
				self.prevN == self.n):
				
				# save coefs in hash table for structure factor
				self.coefs['[%i,%i,%i]'%(g1,g2,g3)] = [re,im]
				
				# save data to list for future
				self.data.append([self.n,		# band number
					self.k1,self.k2,self.k3,	# k-pt
					g1,g2,g3,					# G-pt
					re,im,						# PW coefficient
					self.energy])				# energy
					
			# else, append current dataset to master set
			else:
				# for the large datasets
				self.big_data.append(array(self.data, dtype='f'))
				self.big_coefs.append(self.coefs)
					
				# first entry in new dataset
				self.coefs = {}
				self.coefs['[%i,%i,%i]'%(g1,g2,g3)] = [re,im]
				
				self.data = []
				self.data.append([self.n,
					self.k1,self.k2,self.k3,
					g1,g2,g3,
					re,im,
					self.energy])
					
				self.prevKx = self.k1
				self.prevKy = self.k2
				self.prevKz = self.k3
				self.prevN = self.n
			
			if self.write_out == True:	
				# write data to outfile for mathematica import
				self.out.write("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n" % (self.n,
					self.k1,self.k2,self.k3,
					g1,g2,g3,
					re,im,
					self.energy))
			
			nps -= 1
		
		
	# W/in a kpt, read all band data
	def __readBand(self):
	
		nbds = self.nbands
		while nbds >= 1:
		
			nline = self.gcoeff.readline()
			(n, npws) = nline.split()
			self.n = float(n)
			self.npws = int(npws)
			
			# first iteration, declare the 'previous' values
			if self.start:
				self.prevKx = self.k1
				self.prevKy = self.k2
				self.prevKz = self.k3
				self.prevN = self.n
				self.start = False
				
			# check occupancy
			occline = self.gcoeff.readline()
			self.energy = float(occline.split()[1])
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
		
		if self.readH == False:
			self.readHeader()
			self.readH = True
			
		# open outfile
		if self.write_out == True:
			self.out = open(OUTFILE,"w")
			
		print "Cleaning data..."
		nks = self.nkpts
		while nks >= 1:
		
			k = self.gcoeff.readline()
			(k1, k2, k3) = k.split()
			self.k1 = float(k1)
			self.k2 = float(k2)
			self.k3 = float(k3)
				
			if self.cartK == True:
				self.__getK() # form the actual k-vector from B
				
			self.__readBand()
			nks -= 1
		
		# store data as numpy array
		self.big_data = array(self.big_data)
		
		# close outfile
		if self.write_out == True:
			self.out.close()
		
		self.read = True
			
	# form the structure factor, output to file
	def writeSF(self):
		
		if self.read == False:
			self.readFile()
			self.read = True
		sffile = open(SFFILE,"w")
		print "Writing structure factors..."
		
		for h in range(-5,6,1):
			for k in range(-5,6,1):
				for l in range(-5,6,1):
					coef = self.__getSF(h,k,l)
					sffile.write("%i, %i, %i, %f, %f\n" %(h,k,l,coef[0],coef[1]))

		sffile.close()
		
	# print out the square of the wavefunctions
	def norm(self):
	
		if self.read == False:
			self.readFile()
			self.read = True
			
		for seg in self.big_data:
			n = sum(seg[:,7]**2 + seg[:,8]**2)
			print "%f {%f, %f, %f}:  %f" % (
				seg[:,0][0], seg[:,1][0], seg[:,2][0], seg[:,3][0], n)
		
		
# use the class			

def main():

	filename = sys.argv[1]
	obj = readGCOEFF(filename)
	obj.write_out = False
	obj.readFile()

	start = time.time()		

	obj.writeSF()
	
	end = time.time()

	print "Time elapsed:  %4f"%(end-start)
	
main()

