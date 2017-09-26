#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include <utility>
#include <iostream>
#include <vector>

using namespace std;

int plot_1photon_pT( string histfile="", string histname="none", bool writeplots = true )
{

  /* Default names for debugging */
  histfile="../data/data/DirectPhotonPP-Run13pp510ERT.root";
  histname="pt_1photon";

  //  gStyle->SetOptStat(0);

  /* Open histogram file */
  TFile *f_in = new TFile( histfile.c_str(), "OPEN" );
  f_in->ls();

  THnSparse* h_pT = (THnSparse*) f_in->Get( histname.c_str() );
  cout << "Entries (all pT bins): " << h_pT->GetEntries() << endl;

  /* set cut for direct photons */
  h_pT->GetAxis(2)->SetRange( 0 , 0 );

  TH1F* h_pT_project = (TH1F*)h_pT->Projection(0);

  h_pT_project->Sumw2();
  h_pT_project->GetXaxis()->SetRange(0,30);

  /* Plot pT spectrum */
  TCanvas *c1 = new TCanvas();
  c1->SetLogy();
  h_pT_project->Draw("");

  c1->Print("plots-directphoton/single_photon_pT.eps");
  c1->Print("plots-directphoton/single_photon_pT.png");

  return 0;
}
