#include "GlobalVars.h"
#include "QueryTree.h"

void draw_ERTEff_Photon(const int iso = 0)
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  QueryTree *qt_ert = new QueryTree(Form("data/ERTEff-%sphoton.root",iso?"iso":""), "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-DC3sigma.root");

  TFile *f_pisa = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/HadronResponse-histo-photon.root");
  THnSparse *hn_1photon = (THnSparse*)f_pisa->Get("hn_1photon");

  int evtype = 2;
  int checkmap = 1;
  int ival = 1;

  TH1 *h_ert_t = (TH1*)f->Get("h_ert_0");
  h_ert_t->Reset();
  for(int part=0; part<3; part++)
  {
    TH1 *h_ert[2];
    for(int ert_trig=0; ert_trig<2; ert_trig++)
    {
      h_ert[ert_trig] = (TH1*)h_ert_t->Clone(Form("h_ert_trig%d",ert_trig));
      for(int isolated=iso; isolated<2; isolated++)
      {
        int ih = part + 3*ert_trig + 3*2*evtype + 3*2*3*checkmap + 3*2*3*2*isolated + 3*2*3*2*2*ival;
        TH1 *h2_tmp = (TH1*)f->Get(Form("h_ert_%d",ih));
        h_ert[ert_trig]->Add(h2_tmp);
        delete h2_tmp;
      }
    }
    h_ert[0]->Add(h_ert[1]);
    qt_ert->Fill(h_ert[1], h_ert[0], part);
    delete h_ert[0];
    delete h_ert[1];
  }

  for(int part=0; part<3; part++)
  {
    TH1 *h_ert[2];
    for(int ert_trig=0; ert_trig<2; ert_trig++)
    {
      hn_1photon->GetAxis(1)->SetRange(secl[part],sech[part]);
      hn_1photon->GetAxis(4)->SetRange(ert_trig+1,2);
      h_ert[ert_trig] = (TH1*)hn_1photon->Projection(0);
      h_ert[ert_trig]->SetName(Form("h_ert_trig%d",ert_trig));
    }
    qt_ert->Fill(h_ert[1], h_ert[0], part+3);
    delete h_ert[0];
    delete h_ert[1];
  }

  mc();
  mcd();
  legi(0, 0.35,0.2,0.55,0.4);

  for(int part=0; part<3; part++)
  {
    TGraphAsymmErrors *gr = qt_ert->GraphAsymm(part);
    gr->SetTitle("ERT_4x4c trigger efficeincy for photon");
    aset(gr, "p_{T} (GeV/c)","Eff", 1.,30., 0.,1.1);
    style(gr, part+20, part+1);
    if(part==0)
      gr->Draw("APE");
    else
      gr->Draw("PE");
    leg0->AddEntry(gr, pname[part], "LPE");

    gr->Fit("pol0", "Q","", 10.,30.);
    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr->FindObject("stats");
    st->SetY1NDC(0.4-part*0.1);
    st->SetY2NDC(0.5-part*0.1);

    TGraphAsymmErrors *gr_sim = qt_ert->GraphAsymm(part+3);
    style(gr_sim, part+24, part+1);
    //gr_sim->Draw("PE");
  }
  leg0->Draw();

  c0->Print(Form("plots/ERTEff-%sphoton.pdf",iso?"iso":""));
  qt_ert->Save();
}
