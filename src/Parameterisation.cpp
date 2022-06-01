#include "Parameterisation.h"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

Parameterisation::Parameterisation( G4double dx, G4double dy, G4double dz ) :
  m_dx( dx ), m_dy( dy ), m_dz( dz )
{
}

void Parameterisation::ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const
{
  G4double Xposition = ( copyNo % 5 ) * m_dx;
  G4double Yposition = ( copyNo % 5 ) * m_dy;
  G4double Zposition = ( copyNo % 5 ) * m_dz;

  if ( copyNo > 4 ) Xposition += 100 * cm;

  G4ThreeVector origin( Xposition, Yposition, Zposition );
  physVol->SetTranslation( origin );
  physVol->SetRotation( 0 );
}
