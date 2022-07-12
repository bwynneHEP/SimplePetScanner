#include "DetectorConstruction.h"
#include "ActionInitialization.h"
#include "DecayTimeFinderAction.h"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "QGSP_BERT_HP.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


int main( int argc, char* argv[] )
{
  // Start interactive session using the command line arguments
  G4UIExecutive* ui = new G4UIExecutive( argc, argv );

  // Parse command-line arguments
  bool useGUI = false;
  int nEvents = 0;
  std::string detectorName = "";
  for ( int argi = 1; argi < argc; ++argi )
  {
    std::string argument = argv[ argi ];

    // Find argument parameter
    std::string nextArgument = "";
    if ( argi + 1 < argc )
    {
      nextArgument = argv[ argi + 1 ];
      if ( nextArgument[0] == '-' ) nextArgument = "";
      else ++argi;
    }

    // Find integer argument
    int nextInteger = 0;
    if ( nextArgument.size() )
    {
      nextInteger = atoi( nextArgument.c_str() );
    }

    // Examine arguments
    if ( argument == "--gui" ) useGUI = true;
    else if ( argument == "-n" )
    {
      if ( nextInteger ) nEvents = nextInteger;
      else
      {
        std::cerr << "Failed to parse number of events" << std::endl;
        return 1;
      }
    }
    else if ( argument == "--detector" )
    {
      if ( nextArgument.size() ) detectorName = nextArgument;
      else
      {
        std::cerr << "Did not find detector name to use" << std::endl;
        return 1;
      }
    }
  }

  // Set up the random number generator
  G4Random::setTheEngine( new CLHEP::RanecuEngine );
  G4Random::setTheSeed( 1234 );

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager();

  // Set up physics processes
  G4VModularPhysicsList* physicsList = new QGSP_BERT_HP();
  physicsList->RegisterPhysics( new G4StepLimiterPhysics() );
  physicsList->RegisterPhysics( new G4RadioactiveDecayPhysics() ); // For the tracers
  runManager->SetUserInitialization( physicsList );

  // Set user action classes
  DecayTimeFinderAction * decayTimeFinder = new DecayTimeFinderAction();
  ActionInitialization * actions = new ActionInitialization( decayTimeFinder );
  runManager->SetUserInitialization( actions );

  // Set up detector
  DetectorConstruction * detector = new DetectorConstruction( decayTimeFinder, detectorName );
  runManager->SetUserInitialization( detector );

  G4VisManager* visManager = nullptr;
  if ( useGUI )
  {
    // Set up display
    visManager = new G4VisExecutive();
    visManager->Initialize();
  }

  // Set up the command line interface
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand( "/run/initialize" );
  if ( useGUI ) UImanager->ApplyCommand( "/control/execute vis.mac" );
  UImanager->ApplyCommand( "/control/execute run.mac" );
  UImanager->ApplyCommand( "/run/beamOn " + std::to_string( nEvents ) ); // even if it's zero, useful to initialise physics
  if ( useGUI ) ui->SessionStart();
  delete ui;

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  //
  if ( useGUI ) delete visManager;
  delete runManager;

  return 0;
}
