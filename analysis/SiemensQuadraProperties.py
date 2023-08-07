import math

import PhysicsConstants as PC

def DetectorRadius():
  return 400 # 82cm radius - 1cm crystal half-depth. Maybe 82 should be inner radius and need to adjust simulation?

def CrystalVolume():
  return 0.32 * 0.32 * 2.0

def CrystalMass( DetectorMaterial ):
  if DetectorMaterial not in PC.densities:
    raise RuntimeError('No density value for this material in PhysicsConstants. Wrong material name?')
  return CrystalVolume() * PC.densities[DetectorMaterial]

def DetectorMass( DetectorMaterial ):
  return CrystalMass(DetectorMaterial) * 243200.0

def DetectorDiscreteLength( Length ):
  nRings = float( math.ceil( Length / 32.0 ) )
  return nRings * 32.0

def DetectorMassLength( Length, DetectorMaterial ):
  nRings = float( math.ceil( Length / 32.0 ) )
  return CrystalMass(DetectorMaterial) * 7600.0 * nRings

def LSOunitsInMass( Mass ):
  LSOmoleGram = (2.0 * 175.0) + 28.1 + (5.0 * 16.0) #Lu 2 si 1 O 5
  NAvogadro = 6.022E23
  return NAvogadro * Mass / LSOmoleGram

def Lu176atomsInMass( Mass ):
  return LSOunitsInMass( Mass ) * 2.0 * 0.026 #2.6% Lu176

def Lu176decaysInMass( Mass ):
  Lu176halfYears = 3.76E10
  Lu176halfSec = Lu176halfYears * 3.154E7
  return Lu176atomsInMass( Mass ) * math.log(2.0) / Lu176halfSec
