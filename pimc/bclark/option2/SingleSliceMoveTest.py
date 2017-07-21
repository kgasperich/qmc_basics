from PIMC import *
#import numpy as np
tau=0.1
lam=0.5
numTimeSlices=5
numParticles=2
Path=PathClass(np.zeros((numTimeSlices,numParticles,3),float),tau,lam)
Path.SetPotential(ZeroFunction)
Path.SetCouplingConstant(0.0)
NumSteps=2000
print(PIMC(NumSteps,Path,[SingleSliceMove]))
