#include "GlobalVars.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void draw_ERTEff_Pion()
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  TGraphAsymmErrors *gr[2];
  int igp[2] = {};
  for(int part=0; part<2; part++)
  {
    gr[part] = new TGraphAsymmErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
    for(int ic=0; ic<2; ic++)
      mc(part*2+ic, 6,5);
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[part][cond]
  TH2 *h2_ert_pion[2][2];

  int bbc10cm = 1;
  int ert_trig[2] = {2, 5};

  TH2 *h2_ert_pion_t = (TH2*)f->Get("h2_ert_pion_0");
  h2_ert_pion_t->Reset();
  for(int part=0; part<2; part++)
  {
    for(int cond=0; cond<2; cond++)
    {
      h2_ert_pion[part][cond] = (TH2*)h2_ert_pion_t->Clone(Form("h2_ert_pion_part%d_cond%d",part,cond));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
      {
        int ih = sector + 8*ert_trig[cond] + 6*8*bbc10cm;
        TH2 *h2_tmp = (TH2*)f->Get(Form("h2_ert_pion_%d",ih));
        h2_ert_pion[part][cond]->Add(h2_tmp);
        delete h2_tmp;
      }
    }
    h2_ert_pion[part][0]->Add(h2_ert_pion[part][1]);
  }

  for(int part=0; part<2; part++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      TH1 *h_minv;

      double nt, ent;
      mcd(part*2, ipt+1);
      h_minv = h2_ert_pion[part][0]->ProjectionY("h_minv", ipt+1,ipt+1);
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, nt, ent);
      delete h_minv;

      double np, enp;
      mcd(part*2+1, ipt+1);
      h_minv = h2_ert_pion[part][1]->ProjectionY("h_minv", ipt+1,ipt+1);
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, np, enp);
      delete h_minv;

      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double yy, eyyl, eyyh;
      if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
      {
        eyyl = yy * sqrt( pow(ent/nt,2) + pow(enp/np,2) );
        eyyh = 0.;
      }
      if( yy >= 0. && eyyl >= 0. && eyyl < TMath::Infinity() )
      {
        gr[part]->SetPoint(igp[part], xx, yy);
        gr[part]->SetPointError(igp[part], 0.,0., eyyl,eyyh);
        igp[part]++;
      }
    }

  mc(4);
  mcd(4);

  for(int part=0; part<2; part++)
  {
    gr[part]->SetName(Form("gr_%d",part));
    gr[part]->SetTitle("ERT_4x4c trigger efficeincy for #pi^{0}");
    aset(gr[part], "p_{T} [GeV]","Eff", 1.,30., 0.,1.1);
    style(gr[part], part+20, part+1);
    if(part==0)
      gr[part]->Draw("APE");
    else
      gr[part]->Draw("PE");

    gr[part]->Fit("pol0", "Q","", 10.,20.);
    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr[part]->FindObject("stats");
    st->SetY1NDC(0.3-part*0.2);
    st->SetY2NDC(0.5-part*0.2);

    TGraph *gr_sasha =  new TGraph( Form("data/sasha-trig-part%d.txt",part) );
    gr_sasha->Draw("C");
  }

  legi(0, 0.7,0.2,0.9,0.4);
  leg0->AddEntry(gr[0], "PbSc", "LPE");
  leg0->AddEntry(gr[1], "PbGl", "LPE");
  //leg0->Draw();

  c4->Print("plots/ERTEff-pion.pdf");

  TFile *f_out = new TFile("data/ERTEff-pion.root", "RECREATE");
  for(int part=0; part<2; part++)
  {
    mcw( part*2, Form("part%d-total",part) );
    mcw( part*2+1, Form("part%d-passed",part) );
    gr[part]->Write();
  }
  f_out->Close();
}
