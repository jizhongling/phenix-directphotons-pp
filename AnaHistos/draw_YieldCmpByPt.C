#include "GlobalVars.h"

void draw_YieldCmpByPt()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  TGraph *gr[3];
  Int_t igp[3] = {};
  for(Int_t part=0; part<3; part++)
    gr[part] =  new TGraph(25);

  TFile *f_mine = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  TFile *f_sasha = new TFile("data/Pi0PP-histo.root");

  THnSparse *hn_pion = (THnSparse*)f_mine->Get("hn_pion");
  TAxis *axis_sec = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_cut = hn_pion->GetAxis(3);
  TAxis *axis_type = hn_pion->GetAxis(4);

  for(Int_t part=0; part<3; part++)
    for(Int_t ipt=2; ipt<25; ipt++)
    {
      axis_type->SetRange(3,3);
      axis_cut->SetRange(4,4);
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);

      TH1 *h_minv = hn_pion->Projection(2);
      Double_t npion_mine = h_minv->Integral(113,162);
      delete h_minv;

      TH1 *mchist = (TH1*)f_sasha->Get(Form("mchist_s%d_pt%02d_tp",part,ipt));
      Double_t npion_sasha = mchist->Integral(113,162);
      delete mchist;

      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t ratio = npion_mine / npion_sasha;
      gr[part]->SetPoint(igp[part], xx, ratio);
      igp[part]++;
    }

  mc();
  mcd();
  legi();

  for(Int_t part=0; part<3; part++)
  {
    gr[part]->Set(igp[part]);
    aset(gr[part], "p_{T} [GeV]","#frac{Mine}{Sasha}", 0.,20., 0.49999,0.50001);
    style(gr[part], part+24, part+1);
    if(part == 0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
    leg0->AddEntry(gr[part], pname[part], "P");
  }
  leg0->Draw();

  c0->Print("plots/YieldCmpByPt.pdf");
}
