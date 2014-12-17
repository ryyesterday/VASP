#!/usr/bin/python
# Ryan Valenza
# 2014-12-16
# This script is part of a project to create a series of modules
# to interpret VASP outfiles.  The end goal is to be able to quickly
# create plots and format data for external use.

# Imports
import numpy as np
import math
from scipy import constants

# Globals
HBAR = constants.physical_constants['Planck constant over 2 pi in eV s'][0]
C = constants.physical_constants['speed of light in vacuum'][0]
E0 = 10e3 
LAMBDA = (2*math.pi*HBAR*C/E0)*1e10
GAMMA = LAMBDA/(4*math.pi)
H_MIN = 6
H_MAX = 6
K_MIN = 6
K_MAX = 6
L_MIN = 6
L_MAX = 6

# VASP class whose objects are individual VASP runs
class vasp:

	def __init__(self, dir):
		
		# directory where VASP outfiles are kept
		self.dir = dir
		
		# default filenames, can be overwritten
		self.chgcarf = "CHGCAR_sum"
		self.outfile = "PY_STRFAC.txt"
		
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
		a = np.loadtxt(self.chgcar)
		b = a.reshape(nx*ny*nz) # flatten array		
		self.rho = b.reshape((nx,ny,nz),order = 'F')
		self.SF = np.fft.fftn(self.rho)/(nx*ny*nz) # structure factor

	def writeSF(self):
	
		out = open(self.dir + self.outfile,"w")
		print "Writing SFs for %s..."%self.system
		for h in range(-H_MIN,H_MAX + 1,1):
			for k in range(-K_MIN,K_MAX + 1,1):
				for l in range(-L_MIN,L_MAX + 1,1):
					sf = self.SF[h][k][l]
					r = sf.real
					i = sf.imag
					out.write("%i\t%i\t%i\t%f\t%f\n"%(h,k,l,r,i))
		
		out.close()
		
def main():
	LiF0 = vasp("E:\\Real Space Charge Density\\LiF\\0\\")
	LiF0.readCHGCAR()
	LiF0.writeSF()

	LiF1 = vasp("E:\\Real Space Charge Density\\LiF\\1\\")
	LiF1.readCHGCAR()
	LiF1.writeSF()

	LiF2 = vasp("E:\\Real Space Charge Density\\LiF\\2\\")
	LiF2.readCHGCAR()
	LiF2.writeSF()

	LiF3 = vasp("E:\\Real Space Charge Density\\LiF\\3\\")
	LiF3.readCHGCAR()
	LiF3.writeSF()
	
	LiF4 = vasp("E:\\Real Space Charge Density\\LiF\\4\\")
	LiF4.readCHGCAR()
	LiF4.writeSF()
	
	LiF5 = vasp("E:\\Real Space Charge Density\\LiF\\5\\")
	LiF5.readCHGCAR()
	LiF5.writeSF()
	
	LiF6 = vasp("E:\\Real Space Charge Density\\LiF\\6\\")
	LiF6.readCHGCAR()
	LiF6.writeSF()
	
	LiF7 = vasp("E:\\Real Space Charge Density\\LiF\\7\\")
	LiF7.readCHGCAR()
	LiF7.writeSF()	
	
	LiF8 = vasp("E:\\Real Space Charge Density\\LiF\\8\\")
	LiF8.readCHGCAR()
	LiF8.writeSF()
	
	LiF9 = vasp("E:\\Real Space Charge Density\\LiF\\9\\")
	LiF9.readCHGCAR()
	LiF9.writeSF()
	
	LiF10 = vasp("E:\\Real Space Charge Density\\LiF\\10\\")
	LiF10.readCHGCAR()
	LiF10.writeSF()
	
main()