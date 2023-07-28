import math

import PhysicsConstants as PC

def DetectorRadius():
  return 393

def CrystalVolume():
  return 0.276 * 0.276 * 1.81

def CrystalMass( DetectorMaterial ):
  varName = f'density_{DetectorMaterial}'
  if varName not in PC.__dict__:
    raise RuntimeError(f'{varName} does not exist in PhysicsConstants. Wrong material name?')
  return CrystalVolume() * eval(f'PC.{varName}')

def DetectorVolume():
  return CrystalVolume() * 564480.0

def DetectorMass( DetectorMaterial ):
  return CrystalMass(DetectorMaterial) * 564480.0

def DetectorDiscreteLength( Length ):
  nRings = float( math.ceil( Length / 231.84 ) )
  return nRings * 231.84

def DetectorVolumeLength( Length ):
  nRings = float( math.ceil( Length / 231.84 ) )
  return CrystalVolume() * 70560.0 * nRings

def DetectorMassLength( Length, DetectorMaterial ):
  nRings = float( math.ceil( Length / 231.84 ) )
  return CrystalMass(DetectorMaterial) * 70560.0 * nRings

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
  return 400.0 * Volume
