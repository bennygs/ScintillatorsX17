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
   fCalDetEdepHCID(-1),
   fVBDetEdepHCID(-1),
   fHBDetEdepHCID(-1),
   fCalDetTrackLengthHCID(-1),
   fVBDetTrackLengthHCID(-1),
   fHBDetTrackLengthHCID(-1)
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
                            G4double CalDetEdep,G4double CalDetTrackLength,
                            G4double VBDetEdep, G4double VBDetTrackLength,
                            G4double HBDetEdep, G4double HBDetTrackLength) const
{
  // Print event statistics
  //
  G4cout
     << "   Calorimeter: total energy: "
     << std::setw(7) << G4BestUnit(CalDetEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(CalDetTrackLength, "Length")
     << "   Vertical Bars: total energy: "
     << std::setw(7) << G4BestUnit(VBDetEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(VBDetTrackLength, "Length")
     << "   Horizontal Bars: total energy: "
     << std::setw(7) << G4BestUnit(HBDetEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(HBDetTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
   // Get hist collections IDs
  if ( fCalDetEdepHCID == -1 ) {
    fCalDetEdepHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("CalDet/Edep");
    fVBDetEdepHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("VBDet/Edep");
    fHBDetEdepHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("HBDet/Edep");
    fCalDetTrackLengthHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("CalDet/TrackLength");
    fVBDetTrackLengthHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("VBDet/TrackLength");
    fHBDetTrackLengthHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("HBDet/TrackLength");

  }

  // Get sum values from hits collections
  //
  auto CalDetEdep = GetSum(GetHitsCollection(fCalDetEdepHCID, event));
  auto VBDetEdep = GetSum(GetHitsCollection(fVBDetEdepHCID, event));
  auto HBDetEdep = GetSum(GetHitsCollection(fHBDetEdepHCID, event));

  auto CalDetTrackLength
    = GetSum(GetHitsCollection(fCalDetTrackLengthHCID, event));

  auto VBDetTrackLength
    = GetSum(GetHitsCollection(fVBDetTrackLengthHCID, event));

  auto HBDetTrackLength
    = GetSum(GetHitsCollection(fHBDetTrackLengthHCID, event));

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  //
  analysisManager->FillH1(0, CalDetEdep);
  analysisManager->FillH1(1, VBDetEdep);
  analysisManager->FillH1(2, HBDetEdep);
  analysisManager->FillH1(3, CalDetTrackLength);
  analysisManager->FillH1(4, VBDetTrackLength);
  analysisManager->FillH1(5, HBDetTrackLength);

  // fill ntuple
  //
  analysisManager->FillNtupleDColumn(0, CalDetEdep);
  analysisManager->FillNtupleDColumn(1, VBDetEdep);
  analysisManager->FillNtupleDColumn(2, HBDetEdep);
  analysisManager->FillNtupleDColumn(3, CalDetTrackLength);
  analysisManager->FillNtupleDColumn(4, VBDetTrackLength);
  analysisManager->FillNtupleDColumn(5, HBDetTrackLength);
  analysisManager->AddNtupleRow();

  //print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;
    PrintEventStatistics( CalDetEdep, CalDetTrackLength,
                          VBDetEdep, VBDetTrackLength, HBDetEdep,
                          HBDetTrackLength);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
