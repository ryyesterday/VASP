#!/usr/bin/python
# Ryan Valenza
# 2014-12-16
# This script is part of a project to create a series of modules
# to interpret VASP outfiles.  The end goal is to be able to quickly
# create plots and format data for external use.

# Imports
import numpy as np
import math

# Globals
H_MIN = 6
H_MAX = 6
K_MIN = 6
K_MAX = 6
L_MIN = 6
L_MAX = 6

# VASP class whose objects are individual VASP runs
class vasp:

	# initialize the vasp run object
	def __init__(self, dir):
		
		# directory where VASP outfiles are kept
		self.dir = dir
		
		# default filenames, can be overwritten
		self.chgcarf = "CHGCAR_sum"
		self.outcarf = "OUTCAR"
		self.occf = "OCC.txt"
		self.sff = "PY_STRFAC.txt"
		
	# reads a VASP CHGCAR file and grabs the charge and SF
	def readCHGCAR(self):
		
		# open the file and read the header
		self.chgcar = open(self.dir + self.chgcarf, "r")
		
		self.system = self.chgcar.readline().strip()
		print "Reading CHGCAR for %s..."%self.system
		
		# lattice constant
		self.lc = float(self.chgcar.readline())
		
		# lattice vectors
		a1 = self.chgcar.readline().split()
		a1 = np.array([float(i)*self.lc for i in a1])
		a2 = self.chgcar.readline().split()
		a2 = np.array([float(i)*self.lc for i in a2])
		a3 = self.chgcar.readline().split()
		a3 = np.array([float(i)*self.lc for i in a3])
		self.A = np.array([a1,a2,a3], dtype = 'f8')
		self.findB(self.A) # find reciprocal vectors
		
		# atoms
		self.atoms = self.chgcar.readline().split()
		self.atoms[-1] = self.atoms[-1].strip()
		
		# count per atom 
		natoms = self.chgcar.readline().split()
		self.natoms = [int(i) for i in natoms]
		total = sum(self.natoms)
		
		# direct or cartesian coordinates; basis vectors
		cord = self.chgcar.readline()[0]
		basis = []
		for i in range(total):
			base = self.chgcar.readline().split()
			base = np.array([float(i) for i in base])
			basis.append(base)
		basis = np.array(basis, dtype = 'f8')
		if cord.upper() == 'C':
			self.basis = self.lc * basis
		elif cord.upper() == 'D':
			self.basis = np.dot(self.A, basis.T).T
		self.chgcar.readline()
		
		# grid dimensions
		grid = self.chgcar.readline().split()
		nx = int(grid[0])
		ny = int(grid[1])
		nz = int(grid[2])
		self.grid = [nx,ny,nz]
		
		# the charge density, column-major order	
		extra = (nx*ny*nz % 5) # VASP requires 5 nums per line
		a = np.loadtxt(self.chgcar)
		b = a.reshape(nx*ny*nz + extra) # flatten array	
		if extra != 0:
			extra -= 5
			b = np.delete(b,np.s_[nx*ny*nz:nx*ny*nz+extra])
		self.rho = b.reshape((nx,ny,nz),order = 'F')
		self.SF = np.fft.fftn(self.rho)/(nx*ny*nz) # structure factor

	# reads a VASP OUTCAR file and finds gamma pt occupation
	def getOCC(self):
		
		# open the file and loop through until line match
		self.outcar = open(self.dir + self.outcarf, "r")
		self.occ = open(self.dir + self.occf, "w")
		
		first = True
		gamma = False
		
		print "Finding gamma point occupation..."
		while True:
			line = self.outcar.readline()
			if line == " k-point     1 :       0.0000    0.0000    0.0000\n":
				break
		self.outcar.readline() # skip header
		while True:
			line = self.outcar.readline()
			if line == "\n":
				break
			self.occ.write(line)
				
	# creates an outfile of SF at a given miller index
	def writeSF(self):
	
		out = open(self.dir + self.sff,"w")
		print "Writing SFs for %s..."%self.system
		for h in range(-H_MIN,H_MAX + 1,1):
			for k in range(-K_MIN,K_MAX + 1,1):
				for l in range(-L_MIN,L_MAX + 1,1):
					sf = self.SF[h][k][l]
					r = sf.real
					i = sf.imag
					out.write("%i\t%i\t%i\t%f\t%f\n"%(h,k,l,r,i))
		
		out.close()
		
	# find the reciprocal lattice vectors
	def findB(self, A):
		
		a1 = A[0];
		a2 = A[1];
		a3 = A[2];
		vol = np.dot(a1,np.cross(a2,a3))
		
		b1 = 2*math.pi*(np.cross(a2,a3)/vol)
		b2 = 2*math.pi*(np.cross(a3,a1)/vol)
		b3 = 2*math.pi*(np.cross(a1,a2)/vol)
		self.B = np.array([b1,b2,b3],dtype='f8')
		
	
		
		
# use the class
def main():

	
	C0 = vasp("E:\\WDC\\Graphite\\0\\")
	#C3.readCHGCAR()
	#C3.writeSF()
	C0.getOCC()
	
main()