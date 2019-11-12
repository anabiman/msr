'''
Created on January 27, 2013
Author: Andrew-AM

The SWM polynomial type of coarse grained variables.
'''

import numpy
from scipy.special import legendre
from scipy.linalg import qr
from numpy.linalg import norm, solve

class SpaceWarpingSubsystem:
    """
    A set of CG variables.
    """
    def __init__(self, atoms, kmax):

        self.atoms = atoms
        self.kmax = kmax

        boxboundary = self.atoms.bbox()
        self.box = (boxboundary[1,:] - boxboundary[0,:]) * 1.1 # 110% of the macromolecular box
        self.ref_com = self.atoms.center_of_mass()
        self.basis = self.Construct_Basis(self.atoms.positions - self.ref_com)
        self.ref_coords = self.atoms.positions.copy()

    def ComputeCGInv(self,CG):
        """
        Computes atomic positions from CG positions
        Using the simplest scheme for now
        @param CG: 3*n_cg x 1 array
        @return: a n_atom x 3array
        """
        nCG = cg.shape[0]/3

        x = numpy.dot(self.basis, cg[:nCG])
        y = numpy.dot(self.basis, cg[nCG:2*nCG])
        z = numpy.dot(self.basis, cg[2*nCG:3*nCG])

        return self.ref_com + numpy.array([x,y,z]).T

    def ComputeCG(self, pos):
        """
        Computes CG momenta or positions
        CG = U^t * Mass * var
        var could be atomic positions or velocities
        """
        Utw = self.basis.T * self.atoms.masses

        cg = solve(numpy.dot(Utw, self.basis), numpy.dot(Utw,pos))

        return cg

    def Construct_Basis(self,coords):
        """
        Constructs a matrix of orthonormalized legendre basis functions
        of size Natoms x NCG.
        """
	# normalize coords to [-1,1]
        scaledPos = (coords - coords.mean(axis=0)) / self.box

        # get basis indices
        pindices = self.poly_indexes(self.kmax)

        # grab the masses, and make it a column vector
        basis = numpy.zeros([scaledPos.shape[0], pindices.shape[0]],'f')

        for i in range(pindices.shape[0]):
            k1, k2, k3 = pindices[i,:]
            px = legendre(k1)(scaledPos[:,0])
            py = legendre(k2)(scaledPos[:,1])
            pz = legendre(k3)(scaledPos[:,2])
            basis[:,i] = px * py * pz

        return basis

    def poly_indexes(self, kmax):
        """
        Create 2D array of Legendre polynomial indices with index sum <= psum.

        For example, if psum is 1, the this returns
        [[0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]]
        Note, the sum of each row is less than or equal to 1.
        """
        indices = []

        for n in range(kmax + 1):
            for i in range(n+1):
                for j in range(n+1-i):
                    indices.append([n-i-j, j, i])

        return numpy.array(indices,'i')

    def coarseGrainOP(self):
        Utw = self.basis.T * self.atoms.masses
        Ub = solve(numpy.dot(Utw, self.basis), Utw)

        return Ub.T

    def fineGrainOP(self):
        Utw = self.basis.T * self.atoms.masses
        Ub = solve(numpy.dot(Utw, self.basis), Utw)
        invOP = numpy.dot(Ub, Ub.T)
        return invOP
