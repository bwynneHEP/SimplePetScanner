#include "EventAction.h"

#include "g4csv.hh"

EventAction::EventAction()
{
}

EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction( const G4Event* )
{
}

void EventAction::EndOfEventAction( const G4Event* )
{
  // Finish the ntuple row
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->AddNtupleRow( 0 );
}
