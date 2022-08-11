import math

def CrystalMass():
  crystalCM3 = 0.32 * 0.32 * 2.0
  return crystalCM3 * 7.4

def DetectorMass():
  return CrystalMass() * 243200.0

def DetectorDiscreteLength( Length ):
  nRings = math.ceil( Length / 32.0 )
  return nRings * 32.0

def DetectorMassLength( Length ):
  nRings = math.ceil( Length / 32.0 )
  return CrystalMass() * 7600.0 * nRings

def LSOunitsInMass( Mass ):
  LSOmoleGram = 2.0*175.0 + 28.1 + 5.0*16.0 #Lu 2 si 1 O 5
  NAvogadro = 6.022E23
  return NAvogadro * Mass / LSOmoleGram

def Lu176atomsInMass( Mass ):
  return LSOunitsInMass( Mass ) * 2.0 * 0.026 #2.6% Lu176

def Lu176decaysInMass( Mass ):
  Lu176halfYears = 3.76E10
  Lu176halfSec = Lu176halfYears * 3.154E7
  return Lu176atomsInMass( Mass ) * math.log(2.0) / Lu176halfSec

#https://scipython.com/book2/chapter-6-numpy/examples/simulating-radioactive-decay/
#The probability that a given nucleus will decay in time dt is dt/tau
#tau = tHalf/ln2
