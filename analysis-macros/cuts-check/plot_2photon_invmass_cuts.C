#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include <utility>
#include <iostream>
#include <vector>

using namespace std;

int plot_2photon_invmass_cuts( string histfile="", string histname="none", bool writeplots = true )
{

  /* Default file name */
  histfile="../data/data/DirectPhotonPP-Run13pp510ERT.root";
  histfile="/gpfs/mnt/gpfs02/phenix/spin3/nfeege/taxi_test/keep/DirectPhotonPP-Run13pp510ERT.root";
  histname="inv_mass_2photon";

  //  gStyle->SetOptStat(0);

  /* Define color for each cut */
  int cutcolors[6];
  cutcolors[0] = kBlack;
  cutcolors[1] = kBlack;
  cutcolors[2] = kBlack;
  cutcolors[3] = kBlue;
  cutcolors[4] = kRed;
  cutcolors[5] = kGreen+1;
  cutcolors[6] = kOrange;

  /* Open histogram file */
  TFile *f_in = new TFile( histfile.c_str(), "OPEN" );
  f_in->ls();

  /* Get histogram */
  THnSparse* h_inv_mass_all = (THnSparse*) f_in->Get( histname.c_str() );
  cout << "Entries (all bins): " << h_inv_mass_all->GetEntries() << endl;

  TH1F* h_inv_mass_project[2][17][7]; // [PbSc/PbGl][pTbin][cut]

  TH1F* h_check_counts = (TH1F*)h_inv_mass_all->Projection( 3 );

  /* loop sector */

  /* select sector */
  int jSector=0;
  h_inv_mass_all->GetAxis(2)->SetRange( 1 , 6 ); // PbSc
  //  h_inv_mass_all->GetAxis(2)->SetRange( 7 , 8 ); // PbGl

  /* loop pT bin */
  int bini = 7;
  for ( int jpT = 1; jpT < 17; jpT++ )
    {
      cout << "Plot " << jpT << " => bin " << bini << " from "
	   << h_inv_mass_all->GetAxis(0)->GetBinCenter(bini) - 0.5*h_inv_mass_all->GetAxis(0)->GetBinWidth(bini)
	   << " to "
	   << h_inv_mass_all->GetAxis(0)->GetBinCenter(bini) + 0.5*h_inv_mass_all->GetAxis(0)->GetBinWidth(bini)
	   << endl;

      TString htitle("invMass_PbSc_pTbin");
      htitle+=bini;

      /* select pT range */
      h_inv_mass_all->GetAxis(0)->SetRange( bini , bini );

      /* loop cuts */
      for ( int jCut = 2; jCut < 7; jCut++ )
	{
	  cout << "Cut " << jCut << " from "
	       << h_inv_mass_all->GetAxis(3)->GetBinCenter(jCut) - 0.5*h_inv_mass_all->GetAxis(3)->GetBinWidth(jCut)
	       << " to "
	       << h_inv_mass_all->GetAxis(3)->GetBinCenter(jCut) + 0.5*h_inv_mass_all->GetAxis(3)->GetBinWidth(jCut)
	       << endl;


	  TString hname = htitle;
	  hname+="_cut";
	  hname+=jCut;

	  /* select cut */
	  h_inv_mass_all->GetAxis(3)->SetRange( jCut , jCut );

	  /* Project histogrmas */
	  h_inv_mass_project[jSector][jpT][jCut] = (TH1F*)h_inv_mass_all->Projection( 1 );
	  h_inv_mass_project[jSector][jpT][jCut]->SetName(hname);
	  h_inv_mass_project[jSector][jpT][jCut]->SetTitle(htitle);
	  h_inv_mass_project[jSector][jpT][jCut]->SetLineColor(cutcolors[jCut]);
	}

      bini++;
    }

  /* create legend */
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(h_inv_mass_project[0][1][3],"Emin","l");
  leg->AddEntry(h_inv_mass_project[0][1][4],"Emin && TOF","l");
  leg->AddEntry(h_inv_mass_project[0][1][5],"Emin && shape ","l");
  leg->AddEntry(h_inv_mass_project[0][1][6],"Emin && CV","l");
  leg->AddEntry(h_inv_mass_project[0][1][2],"Emin && TOF && shape && CV","l");

  /* check number of entries */
  TCanvas *c3 = new TCanvas();
  h_check_counts->Draw();


  /* example single pT bin */
  TCanvas *c1 = new TCanvas();
  //c1->SetLogy();
  h_inv_mass_project[0][1][3]->Draw();
  h_inv_mass_project[0][1][4]->Draw("same");
  h_inv_mass_project[0][1][5]->Draw("same");
  h_inv_mass_project[0][1][6]->Draw("same");
  h_inv_mass_project[0][1][2]->Draw("same");

  leg->Draw();

//  h_inv_mass_project[1]->Fit("fit_pi0";)
//
//  c1->Print("photon-plots/two_photon_invariant_mass_pTexample.eps");
//  c1->Print("photon-plots/two_photon_invariant_mass_pTexample.png");
//
//
//  /* Plot multiple pT bins */
//  TCanvas *c2 = new TCanvas();
//  //c2->SetLogy();
//  c2->Divide(4,4);
//  for ( unsigned j = 1; j < 17; j++ )
//    {
//      c2->cd(j);
//      h_inv_mass_project[j]->Draw();
//    }
//
//  c2->Print("photon-plots/two_photon_invariant_mass_pTgrid.eps");
//  c2->Print("photon-plots/two_photon_invariant_mass_pTgrid.png");

  return 0;
}
