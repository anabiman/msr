import sys
import numpy
import MDAnalysis as md
import os
from numpy.linalg import norm, solve
import spacewarping_subsystem as SS

"""
This is a simple script that uses MSR to guide a trajectory of micro (all-atom) states
to new configurations consistent with the imposed coarse-grained (CG) and fine-grained (FG) constraints

The CG method used is space-warping and the code taken from proto_md:
https://github.com/CTCNano/proto_md/tree/master/proto_md
"""

def readCoords(fname, natoms):
        try:
                with open(fname) as f:

                        lines = (line for line in f if not line.strip()[:3].isalpha())
                        r = numpy.loadtxt(lines)
                        pos = numpy.reshape(r, (3, natoms)).T

                return pos

        except Exception:
                raise

def perturb(atoms: md.AtomGroup, noise) -> numpy.ndarray:
        """ Adds noise to the system whilep reserving system COM """
        com = atoms.center_of_mass()
        atoms.positions += (1.0 - 2.0 * numpy.random.random_sample(atoms.positions.shape)) * noise

        return atoms.positions - atoms.center_of_mass() + com

def write(writer, atoms: md.AtomGroup):
        writer.write(atoms)

if __name__ == '__main__':
	# read user input
        _, pdb, tol = sys.argv
        tol = numpy.float(tol)
        clean = True

        # Unperturbed structure
        U = md.Universe(pdb, guess_bonds=True)

        ss = SS.SpaceWarpingSubsystem(U.atoms, kmax=2) 

        bonds = numpy.array(U.bonds.to_indices())
        angles = numpy.array(U.angles.to_indices())
        angles = numpy.array([angles[:,0], angles[:,2]]).T

        indices = numpy.concatenate((bonds, angles), axis=0)
        indices.sort()

        ncons = len(indices)

        leq2 = norm(U.atoms.positions[indices[:,0]] - U.atoms.positions[indices[:,1]], axis=1)**2.0

        numpy.savetxt('Indices.dat', indices, fmt='%d')
        numpy.savetxt('LengthEq.dat', leq2)

        natoms = U.atoms.n_atoms
        nframes = U.trajectory.n_frames
        ocFname = 'output_coords_rec.dat'
        cgFname = 'CG.dat'
        basis = 'basis.dat'
        icFname = 'input_coords_ref.dat'
        invOpFname = 'invOP.dat'

        W = md.Writer('diAlanine.dcd', n_atoms=natoms)

        cgOP = ss.coarseGrainOP()
        numpy.savetxt(basis, cgOP)

        invOP = ss.fineGrainOP()
        numpy.savetxt(invOpFname, invOP)

        positions_orig = U.atoms.positions.copy()
        positions_pert = perturb(U.atoms, noise=2.0)
        U.atoms.positions = positions_pert

        # Write initial frame (perturbed structure)
        write(W, U.atoms)

        # Recover unperturbed atomic positions (for CG construction, etc.)
        U.atoms.positions = positions_orig

        for ts in U.trajectory:

                with open(icFname, 'w') as fp:
                	numpy.savetxt(fp, positions_pert)

                CG = ss.ComputeCG(U.atoms.positions)
                nCG = len(CG)

                numpy.savetxt(cgFname, CG)

                args = f'--natoms {natoms} --ncons {ncons} --ref {icFname} --indices Indices.dat --lengths LengthEq.dat --cg {cgFname} ' +  \
			f'--cgOP {basis} --fgOP {invOpFname} --out {ocFname} --tol {tol} -PC_type jacobi'

                os.system('mpirun -n 1 ../MSR.a ' + args)

                pos = readCoords(ocFname, natoms)
                U.atoms.positions = pos

                write(W, U.atoms)

                os.system('rm {}'.format(cgFname))
                os.system('rm {}'.format(ocFname))
                os.system('rm {}'.format(icFname))

                if ts.frame > 100:
                        break

        W.close()

	# Clean up files if req
        if clean:
                os.system('rm {}'.format(invOpFname))
                os.system('rm {}'.format(basis))
                os.system('rm Indices.dat')
                os.system('rm LengthEq.dat')
                os.system('rm error.dat')
