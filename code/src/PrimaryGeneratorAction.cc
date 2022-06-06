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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4AutoLock.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "Randomize.hh"

#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
/*: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)*/
{
	//defining a GPS (it is used in batch mode)
    fGPS = new G4GeneralParticleSource();


// Definition of a 60Co source

// Gamma 1.173228 MeV
/*  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  G4double size_x = 5.*cm;
  G4double size_y = 5.*cm;
  G4double size_z = 5.*cm;

  G4double x0 = size_x * G4UniformRand();
  G4double y0 = size_y * G4UniformRand();
  G4double z0 = size_z * G4UniformRand();

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fGPS;
  /*delete fParticleGun;*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //  G4double mom_x = 1.-2.*G4UniformRand();
  //  G4double mom_y = 1.-2.*G4UniformRand();
  //  G4double mom_z = 1.-2.*G4UniformRand();

  /*G4double mom_x =1.;
  G4double mom_y =1.;
  G4double mom_z =1.;

// First gamma
G4double E_1 = G4RandGauss::shoot( 1173.228, 2 ) * keV;

// Second gamma
G4double E_2 = G4RandGauss::shoot( 1332.492, 2 ) * keV;


  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName_1;
  G4ParticleDefinition* gamma_1
    = particleTable->FindParticle(particleName_1="gamma");
  fParticleGun->SetParticleDefinition(gamma_1);
  fParticleGun->SetParticleEnergy(E_1);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mom_x,mom_y,mom_z));
  fParticleGun->GeneratePrimaryVertex(anEvent);


    // default particle kinematic
    G4String particleName_2;
    G4ParticleDefinition* gamma_2
      = particleTable->FindParticle(particleName_2="gamma");
    fParticleGun->SetParticleDefinition(gamma_2);
    fParticleGun->SetParticleEnergy(E_2);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mom_x,mom_y,mom_z));
    fParticleGun->GeneratePrimaryVertex(anEvent);*/

fGPS->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
