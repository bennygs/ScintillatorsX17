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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  analysisManager->CreateH1("E_CalDet","Edep in Calorimeter", 20000, 0., 20000*keV);
  analysisManager->CreateH1("E_VBDet","Edep in Vertical Bars", 2000, 0., 2000*keV);
  analysisManager->CreateH1("E_HBDet","Edep in Hotizontal Bars", 2000, 0., 2000*keV);
  analysisManager->CreateH1("L_CalDet","trackL in Calorimeter", 1200, 0., 12000*mm);
  analysisManager->CreateH1("L_VBDet","trackL in Vertical Bars", 1000, 0., 10*mm);
  analysisManager->CreateH1("L_HBDet","trackL in Horizontal Bars", 1000, 0., 10*mm);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("Telescope", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("E_CalDet");
  analysisManager->CreateNtupleDColumn("E_VBDet");
  analysisManager->CreateNtupleDColumn("E_HBDet");
  analysisManager->CreateNtupleDColumn("L_CalDet");
  analysisManager->CreateNtupleDColumn("L_VBDet");
  analysisManager->CreateNtupleDColumn("L_HBDet");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "X17Sim";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl;
    }
    else {
      G4cout << "for the local thread " << G4endl << G4endl;
    }

    G4cout << " E_CalDet : mean = "
       << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;

   G4cout << " E_VBDet : mean = "
      << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;

   G4cout << " E_HBDet : mean = "
      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;

    G4cout << " L_CalDet : mean = "
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;

    G4cout << " L_VBDet : mean = "
      << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl;

    G4cout << " L_HBDet : mean = "
        << G4BestUnit(analysisManager->GetH1(5)->mean(), "Length")
        << " rms = "
        << G4BestUnit(analysisManager->GetH1(5)->rms(),  "Length") << G4endl;


  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
