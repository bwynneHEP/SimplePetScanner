#include "RunAction.h"

#include "g4csv.hh"

RunAction::RunAction() : G4UserRunAction()
{
  // Set number of layers
  G4int const layerNumber = 5;

  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Create histogram of energy distribution
  analysisManager->CreateH1( "LayerTotal", "Layer total energy", layerNumber, 0.5, 0.5 + float( layerNumber ) ); // id 0

  // Render the histograms in a .ps file
  analysisManager->SetH1Plotting( 0, true );

  // Add an ntuple for energy deposits (ntuple id 0)
  analysisManager->CreateNtuple( "Energy", "Deposited energy" );
  analysisManager->CreateNtupleDColumn( "Generated" );
  for ( unsigned int layerIndex = 1; layerIndex <= layerNumber; ++layerIndex )
  {
    std::string columnName = "Layer" + std::to_string( layerIndex );
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
