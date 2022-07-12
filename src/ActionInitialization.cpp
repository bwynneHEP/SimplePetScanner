#include "ActionInitialization.h"
#include "LinearSourceAction.h"
#include "CrystalIntrinsicAction.h"

#include "G4SystemOfUnits.hh"

ActionInitialization::ActionInitialization( DecayTimeFinderAction * decayTimeFinder, std::string sourceName )
  : G4VUserActionInitialization()
  , m_decayTimeFinder( decayTimeFinder )
  , m_sourceName( sourceName )
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::Build() const
{
  if ( m_sourceName.substr( 0, 6 ) == "Linear" ) this->SetUserAction( new LinearSourceAction( -350.0*mm, 350.0*mm, m_sourceName.substr( 6 ) ) );
  else if ( m_sourceName == "Crystal" ) this->SetUserAction( new CrystalIntrinsicAction( -500.0*mm, 500.0*mm, 400.0*mm, 420.0*mm ) );

  this->SetUserAction( m_decayTimeFinder );
}
