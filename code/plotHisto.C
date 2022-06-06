// ROOT macro file for plotting example B4 histograms
//
// Can be run from ROOT session:
// root[0] .x plotHisto.C

{
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // Draw histos filled by Geant4 simulation
  //

  // Open file filled by Geant4 simulation
  TFile f("Silicon_Telescope.root");

  // Create a canvas and divide it into 2x2 pads
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);

  // Draw Eabs histogram in the pad 1
  TH1D* hist1 = (TH1D*)f.Get("E_cry");
  hist1->Draw("HIST");

  // Create a canvas and divide it into 2x2 pads
  TCanvas* c2 = new TCanvas("c2", "", 20, 20, 1000, 1000);

  // Draw Labs histogram in the pad 2
  TH1D* hist2 = (TH1D*)f.Get("L_cry");
  hist2->Draw("HIST");

  // Create a canvas and divide it into 2x2 pads
  TCanvas* c3 = new TCanvas("c3", "", 20, 20, 1000, 1000);

  // Draw Eabs histogram in the pad 1
  TH1D* hist3 = (TH1D*)f.Get("E_Edet");
  hist3->Draw("HIST");

  // Create a canvas and divide it into 2x2 pads
  TCanvas* c4 = new TCanvas("c4", "", 20, 20, 1000, 1000);

  // Draw Labs histogram in the pad 2
  TH1D* hist4 = (TH1D*)f.Get("L_Edet");
  hist4->Draw("HIST");

}
