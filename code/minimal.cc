//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4ScoringManager.hh"
#include "G4UImanager.hh"

#include "DetectorConstruction.hh"
#include "UserActionInitialization.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4Timer.hh"

#include "QGSP_BERT.hh"

#include "Randomize.hh" 
#include <sys/time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
	//Get current time
    G4Timer* theTimer = new G4Timer();
    theTimer->Start();
    
    //Choose the Random Engine 
	G4Random::setTheEngine(new CLHEP::MixMaxRng); //RanecuEngine, MixMaxRng
	time_t systime = time(NULL);
    G4int seed = (long)(systime*G4UniformRand());
    G4Random::setTheSeed(seed);

    //Construct the default run manager
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    int vNumberOfThreads = 1;
    if(argc>2) {
        vNumberOfThreads = atoi(argv[2]);
    }
    if(vNumberOfThreads > 0) {
        runManager->SetNumberOfThreads(vNumberOfThreads);
    }
    G4cout << "### MT MODE ON " << runManager->GetNumberOfThreads() << " ###" << G4endl;
#else
    G4RunManager* runManager = new G4RunManager;
    G4cout << "### MT MODE OFF ###" << G4endl;
#endif

	//Activate UI-command base scorer
	G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
	scManager->SetVerboseLevel(0);
	  
	//Set mandatory initialization classes
	runManager->SetUserInitialization(new QGSP_BERT());
	runManager->SetUserInitialization(new DetectorConstruction);
    runManager->SetUserInitialization(new UserActionInitialization());

	//Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();  

	if (argc!=1) {
		//Batch mode
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);
	}

	else {
		//Define visualization and UI terminal for interactive mode
		G4VisManager* visManager = new G4VisExecutive;
		visManager->Initialize();

		G4UIExecutive* ui = new G4UIExecutive(argc,argv);
		UImanager->ApplyCommand("/control/execute vis.mac");
		ui->SessionStart();
		delete ui;

		delete visManager;
	}

    //Job termination
    delete runManager;
    
    theTimer->Stop();
    G4cout << "Execution terminated" << G4endl;
    G4cout << (*theTimer) << G4endl;
    delete theTimer;
    
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

