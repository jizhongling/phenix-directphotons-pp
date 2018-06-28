#include "GlobalVars.h"

void draw_Correlation()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaPHPythiaDirectPhoton-macros/AnaPHPythia-histo.root");

  THnSparse *hn_corr = (THnSparse*)f->Get("hn_corr");
  TAxis *axis_pt1 = hn_photon->GetAxis(0);
  TAxis *axis_pt2 = hn_photon->GetAxis(1);
  TAxis *axis_dphi = hn_photon->GetAxis(2);
  TAxis *axis_type = hn_photon->GetAxis(3);

  mc(0, 4,4);

  int ipad = 1;
  for(int ipt1=10; ipt1<23; ipt1+=4)
    for(int ipt2=4; ipt2<17; ipt2+=4)
    {
      mcd(0, ipad++);
      for(int itype=0; itype<3; itype++)
      {
        axis_pt1->SetRange(ipt1+1,ipt1+4);
        axis_pt2->SetRange(ipt2+1,ipt2+4);
        axis_type->SetRange(itype+1,itype+1);

        TH1 *h_dphi = hn_corr->Projection(2);
        h_dphi->SetTitle( Form("p_{T}: %.1f-%.1f #otimes %.1f-%.1f GeV",pTbin[ipt1],pTbin[ipt1+4],pTbin[ipt2],pTbin[ipt2+4]) );
        aset(h_dphi);
        style(h_dphi, 20+itype, 1+itype);
        if(itype==0)
          h_dphi->DrawCopy("P0");
        else
          h_dphi->DrawCopy("P0 SAME");
        
        delete h_dphi;
      }
    }

  c0->Print("plots/Correlation.pdf");
}
