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

void draw_Inv_Mass()
{
  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/8511/data/total.root");
  THnSparse *hn_inv_mass_2photon = (THnSparse*)f->Get("inv_mass_2photon");

  TAxis *axis0 = (TAxis*)hn_inv_mass_2photon->GetAxis(0);
  //TAxis *axis1 = (TAxis*)hn_inv_mass_2photon->GetAxis(1);
  TAxis *axis2 = (TAxis*)hn_inv_mass_2photon->GetAxis(2);
  TAxis *axis3 = (TAxis*)hn_inv_mass_2photon->GetAxis(3);

  Int_t Last0 = axis0->GetLast();
  //Int_t Last1 = axis1->GetLast();
  Int_t Last2 = axis2->GetLast();
  Int_t Last3 = axis3->GetLast();

  Int_t nsec = (Last2<8) ? Last2 : 8;
  for(Int_t isec=0; isec<nsec; isec++)
  {
    TCanvas *c = new TCanvas("c", "Canvas", 1200, 1200);
    gStyle->SetOptStat(0);
    c->Divide(4,4);

    axis2->SetRange(isec+1,isec+1);
    Int_t npt = (Last0<15) ? Last0 : 15;
    for(Int_t ipt=1; ipt<=npt; ipt++)
    {
      c->cd(ipt);
      axis0->SetRange(ipt+1,ipt+1);
      Double_t pt_low = axis0->GetBinLowEdge(ipt+1);
      Double_t pt_high = axis0->GetBinLowEdge(ipt+2);
      for(Int_t icon=2; icon<=Last3; icon++)
      {
        axis3->SetRange(icon,icon);
        TH1D *hnp_inv_mass_2photon = (TH1D*)hn_inv_mass_2photon->Projection(1);
        hnp_inv_mass_2photon->SetLineColor(icon);
        if(icon==2)
        {
          char title[100];
          sprintf(title, "Sector_%d_pT_%4.2f_%4.2f", isec, pt_low, pt_high);
          hnp_inv_mass_2photon->SetTitle(title);
          hnp_inv_mass_2photon->DrawCopy();
        }
        else
        {
          hnp_inv_mass_2photon->DrawCopy("SAME");
        }
        hnp_inv_mass_2photon->Delete();
      }
    }

    c->cd(npt+1);
    const char *cond[] = {"Direct Photon", "Photon", "E_{min}", "ToF", "Shape", "#theta_{CV}"};

    Double_t y = 0.8;
    for(Int_t icon=2; icon<=Last3; icon++)
    {

      TLine *line = new TLine(0.2, y, 0.5, y);
      line->SetLineColor(icon);
      line->Draw();

      TLatex *t = new TLatex();
      t->SetTextFont(22);
      t->SetTextAlign(12);
      t->SetNDC();
      t->DrawLatex(0.6, y, cond[icon-1]);

      y -= 0.1;
    }

    char buf[100];
    sprintf(buf, "Inv_Mass-%d.pdf", isec);
    c->Print(buf);
    delete c;
  }
}
