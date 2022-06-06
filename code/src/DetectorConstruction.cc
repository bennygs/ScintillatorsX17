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

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Paraboloid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4AutoDelete.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SDChargedFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true)
  {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VPhysicalVolume* DetectorConstruction::Construct(){

    //** World **//
    G4Material* Galactic = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

    //G4double worldSize = 1. * CLHEP::meter;
    G4double worldSize = 1. * m;

    G4Box* worldSolid = new G4Box("World",
                                  worldSize/2.,
                                  worldSize/2.,
                                  worldSize/2.);

    G4LogicalVolume* worldLogic = new G4LogicalVolume(worldSolid,
                                                      Galactic,
                                                      "World");

    G4PVPlacement* worldPhysical = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     worldLogic,
                                                     "World",
                                                     0,
                                                     false,
                                                     0);


    //** TODO: Insert the material definition here **//

G4Material* Si_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
G4Material* EJ200_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");



    //Print all the materials defined
  	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    G4cout << "###############################" << G4endl;


    //***********//


    //** TODO: Insert the detector geometry here **//

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&  Detector Sizes  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

G4int N_Clovers    =   5;
G4int N_Telescopes =   4;
G4int N_HorBars    =   10;
G4int N_VerBars    =   10;

//  Calorimeter
G4double Cal_x = 10.0*cm;
G4double Cal_y = 5.0*cm;
G4double Cal_z = 5.0*cm;

//***************************** Bars *******************************************

G4double bar_Thickness = 0.2*cm;

//  Horizontal Bars

G4double HorBar_x = bar_Thickness;
G4double HorBar_y = Cal_y;
G4double HorBar_z = Cal_z/N_VerBars;

//  Vertical Bars
G4double VerBar_x = bar_Thickness;
G4double VerBar_y = Cal_y/N_HorBars;
G4double VerBar_z = Cal_z;



//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&  Detector Elements Construction  &&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//==============================================================================


//*************************** Calorimeter **************************************

  auto CalorimeterS
        = new G4Box("Calorimeter",
          Cal_x/2, Cal_y/2, Cal_z/2);

  auto CalorimeterLV
        = new G4LogicalVolume(
          CalorimeterS,                    // its solid
          EJ200_mat,       // its material
          "Cal_LV");                  // its name

//==============================================================================

//***************************** Bars *******************************************

//------------------------------------------------------------------------------
// Horizontal bars

auto HorBarS
    = new G4Box("HorBar",           //its name
    HorBar_x/2, HorBar_y/2, HorBar_z/2); //its size

auto HorBarLV
    = new G4LogicalVolume(
      HorBarS,                      // its solid
      EJ200_mat,                    // its material
      "HorBar_LV");                    // its name

//------------------------------------------------------------------------------

// Vertical bars

auto VerBarS
        = new G4Box("HorBar",           //its name
          VerBar_x/2, VerBar_y/2, VerBar_z/2); //its size

auto VerBarLV
        = new G4LogicalVolume(
          VerBarS,                    // its solid
          EJ200_mat,       // its material
          "VerBar_LV");                  // its name

//==============================================================================

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&  Clover assembly  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


//==============================================================================
//Assemblying the Clover
//==============================================================================


G4AssemblyVolume* assemblyCloverDetectors = new G4AssemblyVolume();

G4double detSou_Dis = 10.*cm;		//Source-Detector Distance

// Translation and Rotation of a clover inside the assembly

   G4ThreeVector T_Cal_0;
   G4ThreeVector T_Cal;
   G4RotationMatrix* R_Cal = new G4RotationMatrix;

   G4ThreeVector T_HorBar_0;
   G4ThreeVector T_HorBar;
   G4RotationMatrix* R_HorBar = new G4RotationMatrix;


   G4ThreeVector T_VerBar_0;
   G4ThreeVector T_VerBar;
   G4RotationMatrix* R_VerBar = new G4RotationMatrix;

 // Rotation of the assembly inside the world
    G4ThreeVector Tm;
    G4RotationMatrix* Rm = new G4RotationMatrix;

// Fill the assembly by the Clovers

    // Calorimeter
    T_Cal_0.setX( HorBar_x + VerBar_x + Cal_x/2. );
    T_Cal_0.setY( 0. );
    T_Cal_0.setZ( 0. );

    R_Cal->rotateX(0.*deg);
    R_Cal->rotateY(0.*deg);
    R_Cal->rotateY(0.*deg);

    //Horizontal bars
    T_HorBar_0.setX(HorBar_x/2. );
    T_HorBar_0.setY( 0. );
    T_HorBar_0.setZ( -( Cal_z - HorBar_z )/2. );

    //Vertical bars
    T_VerBar_0.setX(HorBar_x + VerBar_x/2. );
    T_VerBar_0.setY( -( Cal_y - VerBar_y )/2. );
    T_VerBar_0.setZ( 0. );


    G4double factor = 0.5 * (sqrt(N_Telescopes) - 1.) ;
    G4double theta[] {0., 72., 144., 216., 288.};

    for (int k = 0 ; k < sqrt(N_Telescopes) ; k++ ){

    	for(int j = 0; j < sqrt(N_Telescopes) ; j++ ){


       		T_Cal.setX(T_Cal_0.getX());
       		T_Cal.setY(T_Cal_0.getY() + ((k-factor)*Cal_y) );
       		T_Cal.setZ(T_Cal_0.getZ() + ((j-factor)*Cal_z));

       		assemblyCloverDetectors->AddPlacedVolume(CalorimeterLV, T_Cal, R_Cal);

       		for(int i=0 ; i < N_HorBars ; i++){

       		T_HorBar.setX(T_HorBar_0.getX());
       		T_HorBar.setY(T_HorBar_0.getY() + ((k-factor)*Cal_y));
       		T_HorBar.setZ(T_HorBar_0.getZ()+ i*HorBar_z +  ((j-factor)*Cal_z) );

       		R_HorBar->rotateX(0.*deg);
       		R_HorBar->rotateY(0.*deg);
       		R_HorBar->rotateZ(0.*deg);

       		assemblyCloverDetectors->AddPlacedVolume(HorBarLV, T_HorBar, R_HorBar);
       		}

       		for(int i=0 ; i < N_VerBars ; i++){
       		T_VerBar.setX(T_VerBar_0.getX());
       		T_VerBar.setY(T_VerBar_0.getY()+ i*VerBar_y + ((k-factor)*Cal_y));
       		T_VerBar.setZ(T_VerBar_0.getZ() +  ((j-factor)*Cal_z));

       		R_VerBar->rotateX(0.*deg);
       		R_VerBar->rotateY(0.*deg);
       		R_VerBar->rotateZ(0.*deg);

       		assemblyCloverDetectors->AddPlacedVolume(VerBarLV, T_VerBar, R_VerBar);
       		}

    	}
    }


        // Now instantiate the layers

       for( int i = 0; i < N_Clovers ; i++ ) {
         // Translation of the assembly inside the world

         if(i == 0){
           Tm.setX(detSou_Dis*cos(theta[i]*deg));
           Tm.setY(detSou_Dis*sin(theta[i]*deg));
           Tm.setZ(0.*cm);

           Rm->rotateX(0.*deg);
           Rm->rotateY(0.*deg);
           Rm->rotateZ(theta[i]*deg);

         }else{
           Tm.setX(detSou_Dis*cos(theta[i]*deg));
           Tm.setY(detSou_Dis*sin(theta[i]*deg));
           Tm.setZ(0.*cm);

           Rm->rotateX(0.*deg);
           Rm->rotateY(0.*deg);
           Rm->rotateZ((theta[i] - theta[i-1])*deg);
         }





         assemblyCloverDetectors->MakeImprint( worldLogic, Tm, Rm );
       }

//555555555555555555555555555555555555555555555555555555555555555555555555555555

    G4double dE_t = 0.011 * mm;
    G4double dE_r = 1.0 * cm;

    G4Tubs* crySolid = new G4Tubs("cry",
      0.*cm,   // Inner radius
      dE_r,   // Outer radius
      dE_t,     // Half length in z
      0.*rad,   // Starting phi angle in radians
      2.*pi*rad);  // Angle of the segment in radians


    G4LogicalVolume* cryLogic = new G4LogicalVolume(crySolid,   //its solid
                                                    Si_mat,     //its material
                                                    "cryLV");     //its name



  // Distance from the center

  G4double r_0 = 11.*cm;

  //define a rotation matrix
  G4double tetha  =   165*deg;
  //G4double phi    =   pi/4.*rad;

  G4RotationMatrix* myRotationMatrix = new G4RotationMatrix();
	myRotationMatrix ->  rotateX(0.);
  myRotationMatrix ->  rotateY(180.*deg - tetha);
  myRotationMatrix ->  rotateZ(0.);


  // define a translation vector
  /*G4ThreeVector posCry = G4ThreeVector( r_0*sin(tetha)*cos(phi),
                                        r_0*sin(tetha)*sin(phi),
                                        r_0*cos(phi) );*/

G4ThreeVector posCry = G4ThreeVector( r_0*sin(tetha),
                                      0.,
                                      r_0*cos(tetha));

  new G4PVPlacement(  myRotationMatrix,             //no rotation
                      posCry,         //at (0,0,0)
                      cryLogic,                //its logical volume
                      "cry",              //its name
                      worldLogic,              //its mother  volume
                      false,                   //no boolean operation
                      0);                       //copy number



G4double E_t = 0.300 * mm;
G4double E_r = 1.0 * cm;

G4Tubs* E_det_Solid = new G4Tubs("E_det",
                                 0.*cm,   // Inner radius
                                 E_r,   // Outer radius
                                 E_t,     // Half length in z
                                 0.*rad,   // Starting phi angle in radians
                                 2.*pi*rad);  // Angle of the segment in radians


G4LogicalVolume* E_det_Logic = new G4LogicalVolume(E_det_Solid,   //its solid
                                                   Si_mat,     //its material
                                                   "E_det_LV");     //its name
G4double dE_E_dis = 1.*mm;

G4ThreeVector posEdet = G4ThreeVector( (r_0 + dE_E_dis)*sin(tetha),
                                      0.,
                                      (r_0 + dE_E_dis)*cos(tetha));

new G4PVPlacement(myRotationMatrix,     //no rotation
                  posEdet,              //at (0,0,0)
                  E_det_Logic,          //its logical volume
                  "E_det",              //its name
                  worldLogic,           //its mother  volume
                  false,                //no boolean operation
                  0);                   //copy number

    //***********//

    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField(){

    //** TODO: Insert the sensitive detectors (SD) here **//
    //
    // Scorers
    //

    // declare Absorber as a MultiFunctionalDetector scorer
    //
    G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

    auto cryDetector = new G4MultiFunctionalDetector("Crystal");
    G4SDManager::GetSDMpointer()->AddNewDetector(cryDetector);

    G4VPrimitiveScorer* primitive;
    primitive = new G4PSEnergyDeposit("Edep");
    cryDetector->RegisterPrimitive(primitive);

    primitive = new G4PSTrackLength("TrackLength");
    auto charged = new G4SDChargedFilter("chargedFilter");
    primitive ->SetFilter(charged);
    cryDetector->RegisterPrimitive(primitive);

    SetSensitiveDetector("HorBar_LV",cryDetector);

    //***********//
    // Definition of the Energy detector as a sensitive detector

    auto E_Detector = new G4MultiFunctionalDetector("Edet");
    G4SDManager::GetSDMpointer()->AddNewDetector(E_Detector);

    primitive = new G4PSEnergyDeposit("Edep");
    E_Detector->RegisterPrimitive(primitive);

    primitive = new G4PSTrackLength("TrackLength");
    primitive ->SetFilter(charged);
    E_Detector->RegisterPrimitive(primitive);

    SetSensitiveDetector("VerBar_LV",E_Detector);
    //***********//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
