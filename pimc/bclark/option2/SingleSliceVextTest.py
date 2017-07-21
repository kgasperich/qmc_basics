from PIMC import *
#import numpy as np

tau=0.5
lam=0.5
numTimeSlices=5
numParticles=2
Path=PathClass(numpy.zeros((numTimeSlices,numParticles,3),float),tau,lam)
Path.SetPotential(HarmonicOscillator)
Path.SetCouplingConstant(0.0)
numSteps=5000
#numSteps=100000
moveList = [SingleSliceMove]
print(PIMC(numSteps,Path,moveList,"SingleSlicePotential"))
