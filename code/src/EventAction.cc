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
/// \file B4dEventAction.cc
/// \brief Implementation of the B4dEventAction class

#include "EventAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
 : G4UserEventAction(),
   fCryEdepHCID(-1),
   fEEdepHCID(-1),
   fCryTrackLengthHCID(-1),
   fETrackLengthHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4THitsMap<G4double>*
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection
    = static_cast<G4THitsMap<G4double>*>(
        event->GetHCofThisEvent()->GetHC(hcID));

  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
  G4double sumValue = 0.;
  for ( auto it : *hitsMap->GetMap() ) {
    // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>
    sumValue += *(it.second);
  }
  return sumValue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(
                            G4double cryEdep, G4double cryTrackLength,
                            G4double EEdep,   G4double ETrackLength) const
{
  // Print event statistics
  //
  G4cout
     << "   Crystal: total energy: "
     << std::setw(7) << G4BestUnit(cryEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(cryTrackLength, "Length")
     << "   E_Detector: total energy: "
     << std::setw(7) << G4BestUnit(EEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(ETrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
   // Get hist collections IDs
  if ( fCryEdepHCID == -1 ) {
    fCryEdepHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("Crystal/Edep");
    fEEdepHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("Edet/Edep");
    fCryTrackLengthHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("Crystal/TrackLength");
    fETrackLengthHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("Edet/TrackLength");

  }

  // Get sum values from hits collections
  //
  auto cryEdep = GetSum(GetHitsCollection(fCryEdepHCID, event));
  auto EEdep = GetSum(GetHitsCollection(fEEdepHCID, event));

  auto cryTrackLength
    = GetSum(GetHitsCollection(fCryTrackLengthHCID, event));

  auto ETrackLength
    = GetSum(GetHitsCollection(fETrackLengthHCID, event));

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  //
  analysisManager->FillH1(0, cryEdep);
  analysisManager->FillH1(1, EEdep);
  analysisManager->FillH1(2, cryTrackLength);
  analysisManager->FillH1(3, ETrackLength);

  // fill ntuple
  //
  analysisManager->FillNtupleDColumn(0, cryEdep);
  analysisManager->FillNtupleDColumn(1, EEdep);
  analysisManager->FillNtupleDColumn(2, cryTrackLength);
  analysisManager->FillNtupleDColumn(3, ETrackLength);
  analysisManager->AddNtupleRow();

  //print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;
    PrintEventStatistics(cryEdep, cryTrackLength, EEdep, ETrackLength);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
