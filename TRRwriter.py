import sys
import numpy as np
import MDAnalysis as md
import os

    
U = md.Universe('0_1_10_11_20.gro')
nFrames =  22
natoms = U.atoms.n_atoms

Writer = md.Writer('trajASR.trr', natoms)    # TRR writer
ts = md.coordinates.base.Timestep(natoms)

ts._unitcell = np.ones(9) * 9374

#os.system('./ASR.a -nc {} -ns {} -c CoordsPert.dat -i Indices.dat -l LengthEq.dat -cg CGs.dat -ncg {} -ref basis.dat -o CoordsTmp.dat'.format(natoms, ncons, nCG))

with open('coordsOut.dat') as f:
	lines = (line for line in f if line[0].isdigit())
        r = np.loadtxt(lines)
	
	for n in range(nFrames):
		rT = r[n*3*natoms:(n+1)*3*natoms]
        	ts._pos[:] = np.reshape(rT, (3, natoms)).T


        	ts.frame = n
        	ts.step = n
        	ts.time = n + 1

        	Writer.write(ts)

	Writer.close()
