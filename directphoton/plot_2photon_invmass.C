#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include <utility>
#include <iostream>
#include <vector>

using namespace std;

int plot_2photon_invmass( string histfile="", string histname="none", bool writeplots = true )
{

  /* Default names for debugging */
  histfile="../data/data/DirectPhotonPP-Run13pp510ERT.root";
  histname="inv_mass_2photon";

  //  gStyle->SetOptStat(0);

  /* Open 2D histogram file */
  TFile *f_in = new TFile( histfile.c_str(), "OPEN" );
  f_in->ls();

  THnSparse* h_inv_mass_allpT = (THnSparse*) f_in->Get( histname.c_str() );
  cout << "Entries (all pT bins): " << h_inv_mass_allpT->GetEntries() << endl;

  TF1* fit_pi0 = new TF1("fit_pi0", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x*x" );
  fit_pi0->SetParameter( 0 , 0.3 );
  fit_pi0->SetParameter( 1 , 0.3 );
  fit_pi0->SetParameter( 2 , 0.3 );
  fit_pi0->SetParameter( 3 , 0.3 );
  fit_pi0->SetParameter( 4 , 0.3 );
  fit_pi0->SetParameter( 5 , 0.3 );
  TF1* fit_pi0_gauss = new TF1("fit_pi0_gauss", "[0]*exp(-0.5*((x-[1])/[2])^2)" );
  fit_pi0_gauss->SetRange(0,1);

  TH1F* h_inv_mass_project[17];

  /* set cut for direct photons */
  h_inv_mass_allpT->GetAxis(3)->SetRange( 0 , 0 );

  bini = 7;
  for ( int ploti = 1; ploti < 17; ploti++ )
    {
      cout << "Plot " << ploti << " => bin " << bini << " from "
	   << h_inv_mass_allpT->GetAxis(0)->GetBinCenter(bini) - 0.5*h_inv_mass_allpT->GetAxis(0)->GetBinWidth(bini)
	   << " to "
	   << h_inv_mass_allpT->GetAxis(0)->GetBinCenter(bini) + 0.5*h_inv_mass_allpT->GetAxis(0)->GetBinWidth(bini)
	   << endl;

      TString hname("invMass_pTbin");
      hname+=bini;

      h_inv_mass_allpT->GetAxis(0)->SetRange( bini , bini );

      h_inv_mass_project[ploti] = (TH1F*)h_inv_mass_allpT->Projection( 1 );
      h_inv_mass_project[ploti]->SetName(hname);
      h_inv_mass_project[ploti]->SetTitle(hname);

      bini++;
    }

  /* Duplicate plots for pT bins and do fit */
  TH1F* h_inv_mass_project_fit[17];
  for ( int ploti = 1; ploti < 17; ploti++ )
    {
      h_inv_mass_project_fit[ploti] = (TH1F*)h_inv_mass_project[ploti]->Clone();

      //  h_inv_mass_project_fit[ploti]->Fit("fit_pi0","","",0.105,0.165);
      h_inv_mass_project_fit[ploti]->Fit("fit_pi0","","",0.05,0.300);

      for ( int p = 0; p < 3; p ++ )
	fit_pi0_gauss->SetParameter(p, fit_pi0->GetParameter(p) );

      cout << "Pi0 peak integral: " << fit_pi0_gauss->Integral(0.0,1.0) << endl;

    }
  /* example single pT bin */
  TCanvas *c1 = new TCanvas();
  //c1->SetLogy();
  h_inv_mass_project[1]->Draw();

  c1->Print("plots-directphoton/two_photon_invariant_mass_pTexample.eps");
  c1->Print("plots-directphoton/two_photon_invariant_mass_pTexample.png");


  /* Plot multiple pT bins */
  TCanvas *c2 = new TCanvas();
  //c2->SetLogy();
  c2->Divide(4,4);
  for ( unsigned j = 1; j < 17; j++ )
    {
      c2->cd(j);
      h_inv_mass_project[j]->Draw();
    }

  c2->Print("plots-directphoton/two_photon_invariant_mass_pTgrid.eps");
  c2->Print("plots-directphoton/two_photon_invariant_mass_pTgrid.png");


  /* Plot multiple pT bins with fit*/
  TCanvas *c3 = new TCanvas();
  //c3->SetLogy();
  c3->Divide(4,4);
  for ( unsigned j = 1; j < 17; j++ )
    {
      c3->cd(j);
      h_inv_mass_project_fit[j]->Draw();
      fit_pi0_gauss->Draw("same");
    }

  c3->Print("plots-directphoton/two_photon_invariant_mass_pTgrid_fit.eps");
  c3->Print("plots-directphoton/two_photon_invariant_mass_pTgrid_fit.png");

  return 0;
}
