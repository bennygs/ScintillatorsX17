// ROOT macro file for plotting example B4 ntuple
//
// Can be run from ROOT session:
// root[0] .x plotNtuple.C

{
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // Draw histos filled by Geant4 simulation
  //

  // Open file filled by Geant4 simulation
  TFile f("crystal.root");

  // Create a canvas for the energy spectrum
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);

  // Get ntuple
  TNtuple* ntuple = (TNtuple*)f.Get("crystal");

  // Draw Eabs histogram in the pad 1
  ntuple->Draw("E_cry");

  // Create a canvas for the tracking
  TCanvas* c2 = new TCanvas("c2", "", 20, 20, 1000, 1000);

  // Draw Labs histogram in the pad 2
  ntuple->Draw("L_cry");


}
