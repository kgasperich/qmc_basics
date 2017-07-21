from PIMC import *
#import numpy as np
numParticles=2
numTimeSlices=5
tau=0.5
lam=0.5
Path=PathClass(ReadArray("data/TestPath.dat"),tau,lam)
Path.SetPotential(ZeroFunction)
Path.SetCouplingConstant(0.0)
print("The value of the kinetic action is ", Path.KineticAction(1,2))
print("The value of the kinetic action is ", Path.KineticActionParticle(1,2,0))
print("The value of the kinetic action is ", Path.KineticActionParticle(1,2,1))
