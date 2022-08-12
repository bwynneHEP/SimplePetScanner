import math

def CrystalVolume():
  return 0.276 * 0.276 * 18.1

def CrystalMass():
  return CrystalVolume() * 7.1

def DetectorVolume():
  return CrystalVolume() * 564480.0

def DetectorMass():
  return CrystalMass() * 564480.0

def DetectorDiscreteLength( Length ):
  nRings = float( math.ceil( Length / 8.0 ) )
  return nRings * 8.0

def DetectorVolumeLength( Length ):
  nRings = float( math.ceil( Length / 8.0 ) )
  return CrystalVolume() * 70560.0 * nRings

def DetectorMassLength( Length ):
  nRings = float( math.ceil( Length / 8.0 ) )
  return CrystalMass() * 70560.0 * nRings

def LSOunitsInMass( Mass ):
  yFraction = 0.2
  LSOmoleGram = ((2.0 - yFraction) * 175.0) + (yFraction * 88.9) + 28.1 + (5.0 * 16.0) #Lu 2-x Y x si 1 O 5
  NAvogadro = 6.022E23
  return NAvogadro * Mass / LSOmoleGram

def Lu176atomsInMass( Mass ):
  yFraction = 0.2
  return LSOunitsInMass( Mass ) * (2.0 - yFraction) * 0.026 #2.6% Lu176

def Lu176decaysInMass( Mass ):
  Lu176halfYears = 3.76E10
  Lu176halfSec = Lu176halfYears * 3.154E7
  return Lu176atomsInMass( Mass ) * math.log(2.0) / Lu176halfSec

# does not give consistent answer with above, might be better
def DecaysInVolume( Volume ):
  return 400.0 * Volume:
