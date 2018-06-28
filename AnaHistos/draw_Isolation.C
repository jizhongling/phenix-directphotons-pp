#include "GlobalVars.h"

void draw_Isolation()
{
  //TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaPHPythiaDirectPhoton-macros/AnaPHPythia-histo.root");
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/Isolation-histo.root");

  THnSparse *hn_photon = (THnSparse*)f->Get("hn_photon");
  TAxis *axis_pt = hn_photon->GetAxis(0);
  TAxis *axis_cone = hn_photon->GetAxis(1);
  TAxis *axis_e = hn_photon->GetAxis(2);
  TAxis *axis_iso = hn_photon->GetAxis(3);
  TAxis *axis_prompt = hn_photon->GetAxis(4);

  mc(0, 2,1);
  mc(1, 2,1);

  int ipad = 1;
  for(int ipt=8; ipt<16; ipt+=4)
  {
    TH1 *h_total[2];  // h_total[iso]
    TH1 *h_prompt[2];  // h_prompt[iso]
    for(int iso=1; iso>=0; iso--)
    {
      axis_pt->SetRange(ipt+1,ipt+4);
      axis_iso->SetRange(iso+1,2);
      //axis_cone->SetRange(5,5);
      axis_e->SetRange(5,5);
      axis_prompt->SetRange(1,2);
      h_total[iso] = hn_photon->Projection(1);
      axis_prompt->SetRange(2,2);
      h_prompt[iso] = hn_photon->Projection(1);
    }

    for(int iso=1; iso>=0; iso--)
    {
      mcd(0, ipad);
      TGraphAsymmErrors *gr_prompt = new TGraphAsymmErrors(h_prompt[iso], h_total[iso]);
      gr_prompt->SetName( Form("gr_prompt_%d",iso) );
      gr_prompt->SetTitle( Form("p_{T}: %.1f-%.1f GeV",pTbin[ipt],pTbin[ipt+4]) );
      aset(gr_prompt, "rcone [rad]","#gamma_{prompt&iso}/#gamma_{iso}");
      style(gr_prompt, 20+iso, 1+iso);
      if(iso==1)
        gr_prompt->Draw("AP");
      else
        gr_prompt->Draw("L");
    }

    mcd(1, ipad);
    TGraphAsymmErrors *gr_iso = new TGraphAsymmErrors(h_prompt[1], h_prompt[0]);
    gr_iso->SetTitle( Form("p_{T}: %.1f-%.1f GeV",pTbin[ipt],pTbin[ipt+4]) );
    aset(gr_iso, "rcone [rad]","#gamma_{prompt&iso}/#gamma_{prompt}");
    style(gr_iso, 21, 2);
    gr_iso->Draw("AP");

    for(int iso=1; iso>=0; iso--)
    {
      delete h_total[iso];
      delete h_prompt[iso];
    }
    ipad++;
  }

  c0->Print("plots/Isolation-PISA-Purity.pdf");
  c1->Print("plots/Isolation-PISA-Efficiency.pdf");
}
