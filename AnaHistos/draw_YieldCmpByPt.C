#include "GlobalVars.h"

void draw_YieldCmpByPt()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TGraph *gr[3];
  int igp[3] = {};
  for(int part=0; part<3; part++)
    gr[part] =  new TGraph(25);

  TFile *f_mine = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");
  TFile *f_sasha = new TFile("data/Pi0PP-histo.root");
  //TFile *f_sasha = new TFile("/phenix/spin/phnxsp01/shura/taxi/Run13pp510ERT/5116/data/nt_merged_ert.root");

  // h2_pion[part]
  TH2 *h2_pion[3];

  int evtype = 2;
  int bbc10cm = 1;
  int prob = 1;
  int ival = 1;

  TH2 *h2_pion_t = (TH2*)f_mine->Get("h2_pion_0");
  h2_pion_t->Reset();
  for(int part=0; part<3; part++)
  {
    h2_pion[part] = (TH2*)h2_pion_t->Clone(Form("h2_pion_%d",part));
    for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
      for(int evenodd=0; evenodd<2; evenodd++)
        for(int pattern=0; pattern<3; pattern++)
          for(int isolated=0; isolated<2; isolated++)
            for(int tof=1; tof<3; tof++)
            {
              int ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isolated + 8*2*3*2*tof + 8*2*3*2*3*prob + 8*2*3*2*3*2*evtype + 8*2*3*2*3*2*4*bbc10cm + 8*2*3*2*3*2*4*2*ival;
              TH2 *h2_tmp = (TH2*)f_mine->Get(Form("h2_pion_%d",ih));
              h2_pion[part]->Add(h2_tmp);
              delete h2_tmp;
            }
  }

  for(int part=0; part<3; part++)
    for(int ipt=2; ipt<25; ipt++)
    {
      TH1 *h_minv = h2_pion[part]->ProjectionY("h_minv", ipt+1,ipt+1);
      double npion_mine = h_minv->Integral(113,162);
      delete h_minv;

      TH1 *mchist = (TH1*)f_sasha->Get(Form("mchist_s%d_pt%02d_tp",part,ipt));
      //TH1 *mchist = (TH1*)f_sasha->Get(Form("mc_s%d_bcc0_pt_%03d_tp",part,5*ipt));
      double npion_sasha = mchist->Integral(113,162);
      delete mchist;

      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double ratio = npion_mine / npion_sasha;
      gr[part]->SetPoint(igp[part], xx, ratio);
      igp[part]++;
    }

  mc();
  mcd();
  legi();

  for(int part=0; part<3; part++)
  {
    gr[part]->Set(igp[part]);
    aset(gr[part], "p_{T} [GeV]","#frac{Mine}{Sasha}", 0.,20., 0.499,0.501);
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
