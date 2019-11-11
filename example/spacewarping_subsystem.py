'''
Created on January 27, 2013
Author: Andrew-AM

The SWM polynomial type of coarse grained variables.
'''

import subsystems
import numpy as np
from scipy.special import legendre
from scipy.linalg import qr
from numpy.linalg import norm, solve

class SpaceWarpingSubsystem(subsystems.SubSystem):
    """
    A set of CG variables.
    """
    def __init__(self, system, pindices, select, freq):
        """
        Create a SWM subsystem.
        @param system: an object (typically the system this subsystem belongs to)
        which has the following attributes:
        box: an array describing the periodic boundary conditions.

        note, the system MAY be None, such as when the config is created, so don't
        access it yet.

        @param pindices: a N*3 array of Legendre polynomial indices.
        @param select: a select string used for Universe.selectAtoms which selects
        the atoms that make up this subsystem.
        """
        self.system = system

        # select statement to get atom group
        self.select = select

        # polynomial indices, N_cg x 3 matrix.
        self.pindices = pindices

        # we need to be called with a valid universe
        self.universe_changed(system.universe)
        
        self.cgStep = 0
        
        # How often should the reference struct be updated
        self.Freq_Update = freq

    def universe_changed(self, universe):
        """
        universe changed, so select our atom group
        """
        self.atoms = universe.select_atoms(self.select)

        # check to see if atoms is valid
        if len(self.atoms) <= 0:
            raise ValueError( \
                "The select statement '{}' is syntactically valid but returned a zero sized AtomGroup".format( \
                self.select))

        # call some stuff on the atomGroup to see if its valid
        self.atoms.bbox()

    def frame(self):
        """
        Returns a 3 tuple of CG variables, each one is a row vector of size n_cg
        """
        cg = self.computeCG(self.atoms.positions)
        cgVel = self.computeCGVel(self.atoms.velocities())
        cgFor = self.computeCGForces(self.atoms.forces)

        cg = np.reshape(cg.T, (cg.shape[0]*cg.shape[1]))
        cgVel = np.reshape(cgVel.T, (cgVel.shape[0]*cgVel.shape[1]))
        cgFor = np.reshape(cgFor.T, (cgFor.shape[0]*cgFor.shape[1]))

        return (cg, cgVel, cgFor)

    def translate(self, dCG):
        """
        translates the atomic positions from a given vectory of dCG positions,
        where dCG is a finite change in the CG velocities times dt,
        and then adds the residuals for better accuracy.

        @param CG: a length N_cg 1D array.
        """
        self.atoms.positions += self.computeCGInv(dCG)

    def minimized(self):
        pass

    def md(self):
    	self.cgStep += 1
    	
    def equilibrated(self):
        """
        this is called just after the structure is equilibriated, this is the starting struct
        for the MD runs, this is to calculate basis.
        """
        if self.cgStep%self.Freq_Update == 0:
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

        x = np.dot(self.basis, cg[:nCG])
        y = np.dot(self.basis, cg[nCG:2*nCG])
        z = np.dot(self.basis, cg[2*nCG:3*nCG])

        return self.ref_com + np.array([x,y,z]).T

    def ComputeCG(self, pos):
        """
        Computes CG momenta or positions
        CG = U^t * Mass * var
        var could be atomic positions or velocities
        """
        Utw = self.basis.T * self.atoms.masses

        cg = solve(np.dot(Utw, self.basis), np.dot(Utw,pos))

        return cg
        
    def ComputeCG_Vel(self, vel):
        """
        Computes CG momenta or positions
        CG = U^t * Mass * var
        var could be atomic positions or velocities
        """
        Utw = self.basis.T * self.atoms.masses	
        vel = solve(np.dot(Utw, self.basis), np.dot(Utw,vel))

        return vel

    def ComputeCG_Forces(self, atomic_forces):
        """
        Computes CG forces = U^t * <f>
        for an ensemble average atomic force <f>
        """
        Utw = self.basis.T * self.atoms.masses

        return solve(np.dot(Utw, self.basis), np.dot(Utw, forces))

    def Construct_Basis(self,coords):
        """
        Constructs a matrix of orthonormalized legendre basis functions
        of size Natoms x NCG.
        """
	# normalize coords to [-1,1]        
        scaledPos = (coords - coords.mean(axis=0)) / self.box

        # grab the masses, and make it a column vector
        basis = np.zeros([scaledPos.shape[0], self.pindices.shape[0]],'f')

        for i in range(self.pindices.shape[0]):
            k1, k2, k3 = self.pindices[i,:]
            px = legendre(k1)(scaledPos[:,0])
            py = legendre(k2)(scaledPos[:,1])
            pz = legendre(k3)(scaledPos[:,2])
            basis[:,i] = px * py * pz

        return basis

def poly_indexes(kmax):
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

    return np.array(indices,'i')

def SpaceWarpingSubsystemFactory(system, **args):
    """
    create a list of LegendreSubsystems.
    @param system: the system that the subsystem belongs to, this may be None
                   when the simulation file is created.
    @param selects: A list of MDAnalysis selection strings, one for each
                    subsystem.
    @param args: a list of length 1 or 2. The first element is kmax, and
                 the second element may be the string "resid unique", which can be
                 thought of as an additional selection string. What it does is
                 generate a subsystem for each residue. So, for example, select
                 can be just "resname not SOL", to strip off the solvent, then
                 if args is [kmax, "resid unique"], an seperate subsystem is
                 created for each residue.
    """
    kmax, freq = 0, 1000

    if 'selects' not in args:
        selects = ['all']
    else:
        selects = args['selects']

    try:
        kmax = args['kmax']
    except:
        raise ValueError("invalid subsystem args")

    if 'freq' in args:
        freq = args['freq']

    # test to see if the generated selects work
    [system.universe.select_atoms(select) for select in selects]

    # create the polynomial indices
    pindices = poly_indexes(kmax)

    # the number of CG variables
    # actually, its sufficient to just say nrows * 3 as the
    # number of columns had better be 3.
    ncg = pindices.shape[0] * pindices.shape[1]

    return (ncg, [SpaceWarpingSubsystem(system, pindices, select, freq) for select in selects])
