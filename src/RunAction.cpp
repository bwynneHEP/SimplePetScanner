#include "RunAction.h"

#include "g4csv.hh"

RunAction::RunAction() : G4UserRunAction()
{
  // Set number of crystals
  G4int const crystalNumber = 200;

  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Add an ntuple for energy deposits (ntuple id 0)
  analysisManager->CreateNtuple( "Energy", "Deposited energy" );
  analysisManager->CreateNtupleDColumn( "Generated" );
  for ( unsigned int crystalIndex = 1; crystalIndex <= crystalNumber; ++crystalIndex )
  {
    std::string columnName = "Layer" + std::to_string( crystalIndex );
    analysisManager->CreateNtupleDColumn( columnName );
  }
  analysisManager->FinishNtuple();
}

RunAction::~RunAction()
{
  // Delete analysis manager
  delete G4AnalysisManager::Instance();
}

void RunAction::BeginOfRunAction( const G4Run* )
{
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  analysisManager->OpenFile( "output.csv" );
}

void RunAction::EndOfRunAction( const G4Run* )
{
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Save output data
  analysisManager->Write();
  analysisManager->CloseFile();
}
