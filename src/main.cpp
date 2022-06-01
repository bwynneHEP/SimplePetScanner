#include "DetectorConstruction.h"
#include "ActionInitialization.h"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


int main( int argc, char* argv[] )
{
  // Start interactive session using the command line arguments
  G4UIExecutive* ui = new G4UIExecutive( argc, argv );

  // Set up the random number generator
  G4Random::setTheEngine( new CLHEP::RanecuEngine );
  G4Random::setTheSeed( 1234 );

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager();

  // Set up detector
  runManager->SetUserInitialization( new DetectorConstruction() );

  // Set up physics processes
  G4VModularPhysicsList* physicsList = new FTFP_BERT();
  physicsList->RegisterPhysics( new G4StepLimiterPhysics() );
  runManager->SetUserInitialization( physicsList );

  // Set user action classes (just the generator really)
  runManager->SetUserInitialization( new ActionInitialization() );

  // Set up display
  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();

  // Set up the command line interface
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand( "/control/execute vis.mac" );
  UImanager->ApplyCommand( "/control/execute run.mac" );
  ui->SessionStart();
  delete ui;

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  //
  delete visManager;
  delete runManager;
}
