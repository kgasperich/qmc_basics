import numpy
import numpy as np
import pylab
import CalcStatistics
from Histogram import *

def HarmonicOscillator(pos,mass=1.0,omega=1.0):
    pot = 0.5 * mass * omega**2 * np.dot(pos,pos)
    return pot
#return the harmonic oscillator
# potential for a single particle at position r1


def ZeroFunction(slice):
   return 0.0

def ReadArray(fname):
    io = open(fname, 'r')
    line = io.readline()
    if (line[0] != '#'):
        print('Error in ReadArray!')
        exit()
    splitline = line.split()[1:]
    dim = []
    for i in splitline:
        dim.append (int(i))
    dimtuple = tuple(dim)

    data = numpy.loadtxt(fname)
    return data.reshape(dimtuple)

class PathClass:
    def __init__(self,beads,tau,lam):
        self.tau=tau
        self.lam=lam
        self.beads=beads.copy()
        self.NumTimeSlices=len(beads)
        self.NumParticles=len(beads[0])
        self.NumDimensions=beads.shape[2]
        print(self)
    def __str__(self):
        rep = "I have setup the path with a temperature of %6.4f and %d particles." % \
            (1.0/(self.tau*self.NumTimeSlices),self.NumParticles)
        return rep
    def beta(self):
        return self.tau*self.NumTimeSlices
    def SetCouplingConstant(self,c):
        self.c=c
    def SetPotential(self,externalPotentialFunction):
        self.VextHelper=externalPotentialFunction
    def Vee(self,R):
        # you will write this
        # using self.c
        return 0.0
    def Vext(self,R):
        return self.VextHelper(R)
    def KineticAction(self,slice1,slice2):
        tot = np.sum(
                np.power(
                    np.linalg.norm(
                        self.beads[slice2]-self.beads[slice1],
                        axis=1),
                    2))
        tot /= (4 * self.lam * self.tau)
        # you will fill this in
        return tot
    def KineticActionParticle(self,slice1,slice2,particle):
        tot = np.power(
                    np.linalg.norm(
                        self.beads[slice2,particle]-self.beads[slice1,particle]),
                    2)
        tot /= (4 * self.lam * self.tau)
        # you will fill this in
        return tot
    def PotentialAction(self,slice1,slice2):
        r1 = np.linalg.norm(self.beads[slice1], axis=1)
        r2 = np.linalg.norm(self.beads[slice2], axis=1)
        tot = 0.0
        tot += np.sum([self.Vext(ri) for ri in r1])
        tot += np.sum([self.Vext(ri) for ri in r2])
        tot /= 2.0
        tot *= self.tau
        # you will fill this in
        return tot
    def PotentialActionTotal(self):
        tot = 0.0
        for slice in self.beads:
            for ri in np.linalg.norm(slice, axis=1):
                tot += self.Vext(ri)
        tot *= self.tau
        # you will fill this in
        return tot
    def RelabelBeads(self):
        slicesToShift=random.randint(0,self.NumTimeSlices-1)
        l=list(range(slicesToShift,len(self.beads)))+list(range(0,slicesToShift))
        self.beads=self.beads[l].copy()
    def KineticEnergy(self):
        # computes kinetic energy
        # dx,dy,dz between consecutive slices
        #dxyz = self.beads - np.roll(self.beads,1,axis=0)
        # |dr| between consecutive slices
        #dr = np.linalg.norm(dxyz,axis=2)
        dr = np.linalg.norm(self.beads - np.roll(self.beads,1,axis=0),axis=2)
        # |dR|
        dR = np.linalg.norm(dr,axis=1)
        # sqrt(sum(dR**2))
        dRtot = np.linalg.norm(dR)

        KE=0.0
        KE += 3*self.NumParticles/(2*self.tau)
#        KE -= np.power(
#                np.linalg.norm(
#                    self.beads - np.roll(self.beads,1,axis=0),
#                axis=(0,1,2)),
#                2)/(4*self.lam*self.NumTimeSlices*self.tau**2)
        KE -= (dRtot**2 / (4*self.lam*self.NumTimeSlices*(self.tau**2)))
#        print(KE)
        return KE
    def PotentialEnergy(self):
        # computes potential energy
        PE=0.0
        for islice in range(self.NumTimeSlices):
            for iptcl in range(self.NumParticles):
                PE += self.Vext(self.beads[islice,iptcl,:])
        return PE/(self.NumTimeSlices+0.0)
    def Energy(self):
        return self.PotentialEnergy()+self.KineticEnergy()
    @staticmethod
    def draw_beads_3d(ax,beads):
        """ draw all beads in 3D
        Inputs:
         ax: matplotlib.Axes3D object
         beads: 3D numpy array of shape (nslice,nptcl,ndim)
        Output:
         ptcls: a list of pairs of plot objects. There is ony entry for each particle. Each entry has two items: line representing the particle and text labeling the particle.
        Effect:
         draw all particles on ax """
  
        nslice,nptcl,ndim = beads.shape
        com = beads.mean(axis=0) # center of mass of each particle, used to label the particles only
  
        ptcls = []
        for iptcl in range(nptcl):
            mypos = beads[:,iptcl,:] # all time slices for particle iptcl
            pos = np.insert(mypos,0,mypos[-1],axis=0) # close beads
  
            line = ax.plot(pos[:,0],pos[:,1],pos[:,2],marker='o') # draw particle
            text = ax.text(com[iptcl,0],com[iptcl,1],com[iptcl,2],'ptcl %d' % iptcl,fontsize=20) # label particle
            ptcls.append( (line,text) )
        return ptcls

from numpy import *
def PIMC(numSteps, Path, moveList, name='test'):
    observableSkip=10
    printSkip=1000
    EnergyTrace=[]
    numAccept = zeros((len(moveList)),int)
    configs=[]
#    PD = PathDump(name+'PathDump.h5')
    for steps in range(0,numSteps):
        for mi in range(0,len(moveList)):
            if (moveList[mi](Path)):
                numAccept[mi] += 1
            if steps % observableSkip==0 and steps>1000:
                EnergyTrace.append(Path.Energy())
                if steps % printSkip == 0:
                  print("step:{:10d}   acc. ratio: {:1.3f}".format(steps,numAccept[mi]/steps))
                  configs.append(Path.beads)
#                PairCorrelationFunction(Path,PairHistogram)
#                CalcDensity(Path,DensityHistogram)
#                PD.dump(Path)
#        for mi in range(0,len(moveList)):
#            print('Accept ratio = %1.3f' % ((numAccept[mi]+0.0)/numSteps))
    print(CalcStatistics.Stats(numpy.array(EnergyTrace)))
    pylab.plot(EnergyTrace)
    pylab.savefig(name+"Energy.png")
    with h5py.File("beads.h5","w") as hf:
        dset = hf.create_dataset("beads",data=configs)
    with h5py.File("energy.h5","w") as hf:
        dset = hf.create_dataset("energy",data=EnergyTrace)
#    numpy.save("test1.npz",EnergyTrace)
#    numpy.savez("test2.npz",EnergyTrace)
#    numpy.savez_compressed("test3.npz",EnergyTrace)
#    PairHistogram.plotMeNorm(name+"PairCorrelation.png")
#    DensityHistogram.plotMe(name+"Density.png")
    pylab.show()

#PIMC.py
def SingleSliceMove(Path,sigma=0.1):
    Path.RelabelBeads()
    #add your things here

    #choose particle to move
    particle = np.random.randint(0,Path.NumParticles)
    #evaluate old action
    s_old_k = Path.KineticActionParticle(0,1,particle)
    s_old_k += Path.KineticActionParticle(2,1,particle)
    s_old_v = Path.PotentialAction(1,1)

    s_old = s_old_k + s_old_v

    #save old positions
    old_slice = Path.beads[1,particle].copy()
#    print(old_slice)

    #move slice
    Path.beads[1,particle] += 2*sigma*(np.random.random(old_slice.shape[0]) - 0.5)
#    print(Path.beads[1,particle])
    #evaluate new action
    s_new_k = Path.KineticActionParticle(0,1,particle)
    s_new_k += Path.KineticActionParticle(2,1,particle)
    s_new_v = Path.PotentialAction(1,1)

    s_new = s_new_k + s_new_v

    #evaluate acceptance probability
    p_accept = np.exp( -(s_new - s_old))
    accepted = True
    
    #accept/reject
    if p_accept < np.random.random():
        Path.beads[1,particle] = old_slice
        accepted = False

    #make sure to remember if you reject the move, to restore the old location of the coordinates
    return accepted #accepted # return true if accepted, otherwise return false
def WholeSliceMove(Path,sigma=0.1):
    Path.RelabelBeads()
    #add your things here

    #evaluate old action
    s_old_k = Path.KineticAction(0,1)
    s_old_k += Path.KineticAction(2,1)
    s_old_v = Path.PotentialAction(1,1)

    s_old = s_old_k + s_old_v

    #save old positions
    old_slice = Path.beads[1].copy()

    #move slice
    Path.beads[1] += 2*sigma*(np.random.random(old_slice.shape) - 0.5)
    
    #evaluate new action
    s_new_k = Path.KineticAction(0,1)
    s_new_k += Path.KineticAction(2,1)
    s_new_v = Path.PotentialAction(1,1)

    s_new = s_new_k + s_new_v

    #evaluate acceptance probability
    p_accept = np.exp( -(s_new - s_old))
    accepted = True
    
    #accept/reject
    if p_accept < np.random.random():
        Path.beads[1] = old_slice
        accepted = False

    #make sure to remember if you reject the move, to restore the old location of the coordinates
    return accepted #accepted # return true if accepted, otherwise return false

def DisplaceMove(Path):
    # First, compute the total potential action for all time slices.
    s_old = Path.PotentialActionTotal()
    # This move won't change the kinetic action, so you don't need to compute it
    # Don't forget the link action from the last slice back to the first!
    # Save a copy of the old path
    savePath = Path.beads.copy()
    # Now, create a random vector for displacement
    #delta = 4.0*numpy.array([random.random()-0.5, random.random()-0.5, random.random()-0.5])
    delta = 4.0*(numpy.random.random(Path.NumDimensions) - 0.5)
    # move all the time slices
    Path.beads += delta[np.newaxis,np.newaxis,:]
    # Compute the new potential action
    s_new = Path.PotentialActionTotal()
    # Accept or reject based on the change in potential action
    p_accept = np.exp( -(s_new - s_old))
    accepted = True
    
    #accept/reject
    if p_accept < np.random.random():
        Path.beads = savePath.copy()
        accepted = False
    # Remember to copy back savePath if you reject
    # Remember to return True if you accept and False if you reject
    return accepted
