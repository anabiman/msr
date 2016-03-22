import sys
import numpy as np
import MDAnalysis as md
import os
from numpy.linalg import norm, solve
import proto_md.subsystems.spacewarping_subsystem as SS

def readCoords(fname, natoms):

        try:
                with open(fname) as f:

                        lines = (line for line in f if line[1].isdigit())
                        r = np.loadtxt(lines)

			print r.shape

                        pos = np.reshape(r, (3, natoms)).T

                return pos

        except:
                print "Unexpected error:", sys.exc_info()[0]
                raise

if __name__ == '__main__':
	# read user input
        _, gro, traj, tpr, tol = sys.argv
	tol = np.float(tol)

        U = md.Universe(gro, traj)

        Top = md.Universe(tpr, gro)

        args = {'kmax':2}
        n, ss = SS.SpaceWarpingSubsystemFactory(U, ['protein'], **args)
        s = ss[0]
        s.universe_changed(U)
        s.equilibrated()

        bonds = np.array(Top.bonds.to_indices())
        angles = np.array(Top.angles.to_indices())
        angles = np.array([angles[:,0], angles[:,2]]).T

        indices = np.concatenate((bonds, angles), axis=0)
        indices.sort()

        ncons = len(indices)

        leq2 = norm(Top.atoms.positions[indices[:,0]] - Top.atoms.positions[indices[:,1]], axis=1)**2.0

        np.savetxt('Indices.dat', indices, fmt='%d')
        np.savetxt('LengthEq.dat', leq2)

        natoms = U.atoms.n_atoms
        nframes = U.trajectory.n_frames
        fname = 'coords.dat'
        cgFname = 'CG.dat'
        basis = 'basis.dat'
        cFname = 'coordsInput.dat'
	invOpFname = 'invOP.dat'

        W = md.Writer('MSR.dcd', n_atoms=natoms)

	Utw = s.basis.T * s.atoms.masses

	Ub = solve(np.dot(Utw, s.basis), Utw)
	np.savetxt(basis, Ub.T)

	invOP = np.dot(Ub, Ub.T)
	np.savetxt(invOpFname, invOP)

        for ts in U.trajectory:

                with open(cFname, 'w') as fp:
                	np.savetxt(fp, s.atoms.positions)

                CG = s.ComputeCG(s.atoms.positions)
                nCG = len(CG)

                np.savetxt(cgFname, CG)

		print 'mpirun -n 1 MSR.a -nc {} -ns {} -c {} -i Indices.dat -l LengthEq.dat -cg {} -ncg {} -ref {} -refTrans {} -o {} -tol {} -PC_type jacobi'.format(natoms, \
                        ncons, cFname, cgFname, nCG, basis, invOpFname, fname, tol)

                os.system('mpirun -n 1 MSR.a -nc {} -ns {} -c {} -i Indices.dat -l LengthEq.dat -cg {} -ncg {} -ref {} -refTrans {} -o {} -tol {} -PC_type jacobi'.format(natoms, \
			ncons, cFname, cgFname, nCG, basis, invOpFname, fname, tol))

                pos = readCoords(fname, natoms)
                U.atoms.set_positions(pos)

                W.write(ts)

                os.system('rm {}'.format(fname))

        W.close()
