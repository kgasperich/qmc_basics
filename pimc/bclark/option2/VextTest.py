from PIMC import *
#import numpy as np
numParticles=2
numTimeSlices=5
tau=0.5
lam=0.5
Path=PathClass(ReadArray("data/TestPath.dat"),tau,lam)
Path.SetPotential(HarmonicOscillator)
Path.SetCouplingConstant(0.0)
print("The value of the external potential is ", Path.Vext(np.array([0.1,0.3,0.1])))
print("The value of the potential action between slice 1 and 2 (with a harmonic external potential is)", Path.PotentialAction(1,2))
