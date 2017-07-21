from PIMC import *
#import numpy as np

tau=0.5
lam=0.5
numTimeSlices=5
numParticles=2
Path=PathClass(numpy.zeros((numTimeSlices,numParticles,3),float),tau,lam)
Path.SetPotential(HarmonicOscillator)
Path.SetCouplingConstant(0.0)
numSteps=50000
#numSteps=500000
#numSteps=5000000
#print(PIMC(numSteps,Path,[SingleSliceMove,DisplaceMove],"SingleSlice"))
print(PIMC(numSteps,Path,[SingleSliceMove,DisplaceMove,WholeSliceMove],"SingleSlice"))
import math
def coth(x):
    return math.cosh(x)/math.sinh(x)
#print("exact",3/2*1.0*coth(numTimeSlices*tau/2))
print("exact",numParticles*3/2*1.0*coth(numTimeSlices*tau/2))
