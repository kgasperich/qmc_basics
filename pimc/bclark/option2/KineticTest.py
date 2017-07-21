from PIMC import *
#import numpy as np

tau=0.5
lam=0.5
numTimeSlices=5
numParticles=2
Path=PathClass(numpy.zeros((numTimeSlices,numParticles,3),float),tau,lam)
Path.SetPotential(ZeroFunction)
Path.SetCouplingConstant(0.0)
#numSteps=50000
numSteps=100000
print(PIMC(numSteps,Path,[SingleSliceMove],"SingleSlice"))
