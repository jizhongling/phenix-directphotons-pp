#include <TFile.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLine.h>

#include <cstdio>

using namespace std;

void draw_Theta_CV()
{
  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/8511/data/total.root");
  THnSparse *hn_inv_mass_2photon_theta_cv = (THnSparse*)f->Get("inv_mass_2photon_theta_cv");

  TAxis *axis0 = (TAxis*)hn_inv_mass_2photon_theta_cv->GetAxis(0);
  TAxis *axis1 = (TAxis*)hn_inv_mass_2photon_theta_cv->GetAxis(1);
  //TAxis *axis2 = (TAxis*)hn_inv_mass_2photon_theta_cv->GetAxis(2);
  TAxis *axis3 = (TAxis*)hn_inv_mass_2photon_theta_cv->GetAxis(3);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2400);
  gStyle->SetOptStat(0);
  c->Divide(4,4);

  axis3->SetRange(1,6);
  //axis3->SetRange(7,8);

  for(Int_t iE=8; iE<38; iE+=2)
  {
    static Int_t ipad = 1;
    static Int_t icol = 1;
    c->cd(ipad);
    icol = 1;

    axis0->SetRange(iE,iE+1);
    Double_t clusterE_low = axis0->GetBinLowEdge(iE);
    Double_t clusterE_high = axis0->GetBinLowEdge(iE+2);

    for(Int_t ith=2; ith<62; ith+=4)
    {
      axis1->SetRange(ith,ith+3);
      TH1D *hnp_inv_mass_2photon_theta_cv = (TH1D*)hn_inv_mass_2photon_theta_cv->Projection(2);
      Int_t nEntries = hnp_inv_mass_2photon_theta_cv->GetEntries();
      hnp_inv_mass_2photon_theta_cv->Scale(1./nEntries);
      hnp_inv_mass_2photon_theta_cv->SetLineColor(icol);
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
      hnp_inv_mass_2photon_theta_cv->Delete();
      icol++;
    }

    ipad++;
  }

  c->cd(16);
  for(Int_t ith=2; ith<62; ith+=4)
  {
    static Int_t icol = 1;
    static Double_t y = 0.9;
    Double_t theta_cv = axis1->GetBinLowEdge(ith);
    char buf[100];
    sprintf(buf, "%5.3f", theta_cv);

    TLine *line = new TLine(0.2, y, 0.4, y);
    line->SetLineColor(icol);
    line->Draw();

    TLatex *t = new TLatex();
    t->SetTextFont(22);
    t->SetTextAlign(12);
    t->SetNDC();
    t->DrawLatex(0.5, y, buf);

    y -= 0.05;
    icol++;
  }

  c->Print("Theta_CV_PbSc.pdf");
}
