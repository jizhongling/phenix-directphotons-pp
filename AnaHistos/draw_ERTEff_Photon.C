#include "GlobalVars.h"
#include "QueryTree.h"

void draw_ERTEff_Photon()
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  QueryTree *qt_ert = new QueryTree("data/ERTEff-photon.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[part][cond]
  TH1 *h_ert[2][2];

  int bbc10cm = 1;
  int ert_trig[2] = {2, 5};

  TH1 *h_ert_t = (TH1*)f->Get("h_ert_0");
  h_ert_t->Reset();
  for(int part=0; part<2; part++)
  {
    for(int cond=0; cond<2; cond++)
    {
      h_ert[part][cond] = (TH1*)h_ert_t->Clone(Form("h_ert_part%d_cond%d",part,cond));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
      {
        int ih = sector + 8*ert_trig[cond] + 6*8*bbc10cm;
        TH1 *h2_tmp = (TH1*)f->Get(Form("h_ert_%d",ih));
        h_ert[part][cond]->Add(h2_tmp);
        delete h2_tmp;
      }
    }
    h_ert[part][0]->Add(h_ert[part][1]);
    qt_ert->Fill(h_ert[part][1], h_ert[part][0], part);
  }

  mc();
  mcd();

  for(int part=0; part<2; part++)
  {
    TGraphAsymmErrors *gr = qt_ert->GraphAsymm(part);
    gr->SetTitle("ERT_4x4c trigger efficeincy for photon");
    aset(gr, "p_{T} [GeV]","Eff", 1.,30., 0.,1.1);
    style(gr, part+20, part+1);
    if(part==0)
      gr->Draw("APE");
    else
      gr->Draw("PE");

    gr->Fit("pol0", "Q","", 10.,30.);
    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr->FindObject("stats");
    st->SetY1NDC(0.3-part*0.2);
    st->SetY2NDC(0.5-part*0.2);
  }

  legi(0, 0.7,0.2,0.9,0.4);
  leg0->AddEntry(gr[0], "PbSc", "LPE");
  leg0->AddEntry(gr[1], "PbGl", "LPE");
  leg0->Draw();

  c0->Print("plots/ERTEff-photon.pdf");
  qt_ert->Save();
}
