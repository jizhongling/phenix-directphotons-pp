#include <TFile.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLegend.h>

#include <cstdio>

using namespace std;

void draw_Theta_CV_1()
{
  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/8511/data/total.root");
  THnSparse *hn_inv_mass_2photon_theta_cv = (THnSparse*)f->Get("inv_mass_2photon_theta_cv");

  TAxis *axis0 = (TAxis*)hn_inv_mass_2photon_theta_cv->GetAxis(0);
  TAxis *axis1 = (TAxis*)hn_inv_mass_2photon_theta_cv->GetAxis(1);
  //TAxis *axis2 = (TAxis*)hn_inv_mass_2photon_theta_cv->GetAxis(2);
  TAxis *axis3 = (TAxis*)hn_inv_mass_2photon_theta_cv->GetAxis(3);

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

  axis3->SetRange(1,6);
  //axis3->SetRange(7,8);

  Int_t iE = 10;
  Int_t icol = 1;
  c->cd();

  axis0->SetRange(iE,iE);
  Double_t clusterE_low = axis0->GetBinLowEdge(iE);
  Double_t clusterE_high = axis0->GetBinLowEdge(iE+1);

  Int_t ith[3] = {1, 10, 30};
  for(Int_t i=0; i<3; i++)
  {
    axis1->SetRange(ith[i],ith[i]);
    TH1D *hnp_inv_mass_2photon_theta_cv = (TH1D*)hn_inv_mass_2photon_theta_cv->Projection(2);
    Int_t nEntries = hnp_inv_mass_2photon_theta_cv->GetEntries();
    hnp_inv_mass_2photon_theta_cv->Scale(1./nEntries);
    hnp_inv_mass_2photon_theta_cv->SetLineColor(icol);

    Double_t theta_cv = axis1->GetBinLowEdge(ith[i]);
    char buf[100];
    sprintf(buf, "#theta_{CV}=%5.3f", theta_cv);
    leg->SetTextSize(0.02);
    leg->AddEntry(hnp_inv_mass_2photon_theta_cv, buf, "l");
    leg->Draw();

    if(icol==1)
    {
      char title[100];
      sprintf(title, "PbSc_ClusterE_%4.2f_%4.2f", clusterE_low, clusterE_high);
      hnp_inv_mass_2photon_theta_cv->SetTitle(title);
      hnp_inv_mass_2photon_theta_cv->DrawCopy();
    }
    else
    {
      hnp_inv_mass_2photon_theta_cv->DrawCopy("SAME");
    }

    //hnp_inv_mass_2photon_theta_cv->Delete();
    icol++;
  }

  c->Print("Theta_CV_PbSc-1.pdf");
}
