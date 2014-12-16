VASP
====
This project's goal is to use VASP to study the real space charge density upon electronic heating by an x-ray pulse.  It is also desirable to obtain the structure factor (the Fourier transform of the charge density) and create diffraction patterns that can be compared with experimental data acquired from synchrotrons.  During the project, I also became interested in understanding plane wave expansions, band structure calculations, Fermi-Dirac statistics, and other things involved in Density Functional Theory.  Here are a collection of my scripts and LaTeX documents for this project.

Currently this project contains the following scripts:
- readCHGCAR.py
- readGCOEFF.py
- readEIGENVAL.py
- readDOSCAR.py


### readCHGCAR.py

This is the current method I'm using to obtain the structure factor from the charge density that VASP outputs in CHGCAR.  Beware that VASP only outputs the valence charge density, so to obtain the core run VASP with the option LAECHG = .TRUE. in INCAR.  This will create the files AECCAR0, containing the core, and AECCAR2, containing the valence.  Sum them to obtain the total density.  This density is read into readCHGCAR.py as a numpy array in column-major order (VASP's default for writing 3D arrays).  A discrete Fourier transform is performed and the result is a 3D array where the index {h,k,l} corresponds to a miller plane.  See the LaTeX document "Structure Factor FFT" for explanations of the theory.  

### readGCOEFF.py

This part of the project was an attempt to read the WAVECAR file from VASP, the thought being that it would be easy to obtain the structure factor from the wavefunctions given in a plane wave expansion (see the LaTeX document "Structure Factor PW").  To ensure that VASP evaluates the wavefunctions throughout the entire Brillouin zone, the symmetry must be turned off by setting SYM = 0 in INCAR.  I used a FORTRAN script called [WaveTrans](http://www.andrew.cmu.edu/user/feenstra/wavetrans/) written by a group at Carnegie Mellon University.  This script reads the WAVECAR file and produces an outfile called GCOEFF.txt which contains the plane wave coefficients indexed by band number, k-point, and G-vector.  My script readGCOEFF.py reads in the GCOEFF.txt file and reconstructs the structure factor - I was also playing around with python OOP structure and implemented a few other useful tools.  The problem with this method is that only the pseudo-wavefunctions for the valence electrons are included in the WAVECAR file.  It is extra work to include the PAW density coefficients and obtain the all-electron density.  

### readEIGENVAL.py

This script reads in VASP's EIGENVAL file and obtains a data set of {Energy, k-point} pairs for every band calculated.  The data is organized in an outfile to be read into some plotting program.  

### readDOSCAR.py

This script reads in VASP's DOSCAR file for a non spin-polarized calculation and plots the DOS and the integrated-DOS using matplotlib.  

# Author:  Ryan Valenza
