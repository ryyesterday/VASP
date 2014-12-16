import numpy

HLINES = 10
CHGCAR = open("CHGCAR_sum","r")
OUTFILE = open("STRFAC.out","w")

# read the header
i = HLINES
while i >= 0:
	CHGCAR.readline()
	i -= 1
	
(nx,ny,nz) = CHGCAR.readline().split()
nx = int(nx)
ny = int(ny)
nz = int(nz)

a = numpy.loadtxt(CHGCAR)
b = a.reshape(a.shape[0]*a.shape[1])
c = b.reshape((nx,ny,nz),order = 'F')
d = numpy.fft.fftn(c)/(nx*ny*nz)

OUTFILE.write("LINE 1\n")
OUTFILE.write("LINE 2\n")
OUTFILE.write("LINE 3\n")

for h in range(-6,7,1):
	for k in range(-6,7,1):
		for l in range(-6,7,1):
			x = d[h][k][l]
			OUTFILE.write("%i\t%i\t%i\t%f\t%f\n"%(h,k,l,x.real,x.imag))
