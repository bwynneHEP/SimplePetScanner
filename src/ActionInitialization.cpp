#include "ActionInitialization.h"
#include "LinearSourceAction.h"
#include "CrystalIntrinsicAction.h"

#include "G4SystemOfUnits.hh"

ActionInitialization::ActionInitialization( DecayTimeFinderAction * decayTimeFinder, std::string sourceName, G4double detectorLength, G4double phantomLength )
  : G4VUserActionInitialization()
  , m_decayTimeFinder( decayTimeFinder )
  , m_sourceName( sourceName )
  , m_detectorLength( detectorLength )
  , m_phantomLength( phantomLength )
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::Build() const
{
  if ( m_sourceName.substr( 0, 6 ) == "Linear" )
  {
    G4double phantomLength = 350.0*mm;
    if ( m_phantomLength >= 0.0 )
    {
      phantomLength = m_phantomLength * mm / 2.0; // half-lengths
    }
    this->SetUserAction( new LinearSourceAction( -phantomLength, phantomLength, m_sourceName.substr( 6 ) ) );
  }
  else if ( m_sourceName == "Siemens" )
  {
    G4double detectorLength = 512.0*mm;
    if ( m_detectorLength >= 0.0 )
    {
      detectorLength = m_detectorLength * mm / 2.0; // half-lengths
    }
    this->SetUserAction( new CrystalIntrinsicAction( -detectorLength, detectorLength, "LSO", 400.0*mm, 420.0*mm ) );
  }
  else
  {
    std::cerr << "Unrecognised source name: " << m_sourceName << std::endl;
    exit(1);
  }

  this->SetUserAction( m_decayTimeFinder );
}
