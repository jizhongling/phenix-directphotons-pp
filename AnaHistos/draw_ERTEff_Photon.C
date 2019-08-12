#include "GlobalVars.h"
#include "QueryTree.h"

void draw_ERTEff_Photon()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const int partl[2] = {0, 2};
  const int parth[2] = {1, 2};

  QueryTree *qt_ert = new QueryTree("data/ERTEff-photon.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[part][cond]
  TH1 *h_ert[2][2];

  int bbc10cm = 1;
  int evtype = 2;
  int ival = 1;

  TH1 *h_ert_t = (TH1*)f->Get("h_ert_0");
  h_ert_t->Reset();
  for(int part=0; part<2; part++)
  {
    for(int ert_trig=0; ert_trig<2; ert_trig++)
    {
      h_ert[part][ert_trig] = (TH1*)h_ert_t->Clone(Form("h_ert_part%d_trig%d",part,ert_trig));
      for(int ipart=partl[part]; ipart<=parth[part]; ipart++)
        for(int isolated=0; isolated<2; isolated++)
        {
          int ih = ipart + 3*ert_trig + 3*2*evtype + 3*2*3*bbc10cm + 3*2*3*2*isolated + 3*2*3*2*2*ival;
          TH1 *h2_tmp = (TH1*)f->Get(Form("h_ert_%d",ih));
          h_ert[part][ert_trig]->Add(h2_tmp);
          delete h2_tmp;
        }
    }
    h_ert[part][0]->Add(h_ert[part][1]);
    qt_ert->Fill(h_ert[part][1], h_ert[part][0], part);
  }

  mc();
  mcd();
  legi(0, 0.7,0.2,0.9,0.4);

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
    leg0->AddEntry(gr, pname[part], "LPE");

    gr->Fit("pol0", "Q","", 10.,30.);
    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr->FindObject("stats");
    st->SetY1NDC(0.3-part*0.2);
    st->SetY2NDC(0.5-part*0.2);
  }
  leg0->Draw();

  c0->Print("plots/ERTEff-photon.pdf");
  qt_ert->Save();
}
