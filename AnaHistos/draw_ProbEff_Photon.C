#include "GlobalVars.h"

void draw_ProbEff_Photon()
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  TGraphAsymmErrors *gr[2];
  for(int part=0; part<2; part++)
  {
    gr[part] = new TGraphAsymmErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  //TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/PhotonEff-histo.root");

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  TAxis *axis_sec = hn_1photon->GetAxis(0);
  TAxis *axis_pt = hn_1photon->GetAxis(1);
  TAxis *axis_pattern = hn_1photon->GetAxis(2);
  TAxis *axis_cut = hn_1photon->GetAxis(3);
  TAxis *axis_type = hn_1photon->GetAxis(4);

  for(int part=0; part<2; part++)
  {
    axis_type->SetRange(3,3);
    axis_sec->SetRange(secl[part],sech[part]);
    axis_cut->SetRange(1,1);
    TH1 *h_total = hn_1photon->Projection(1);
    axis_cut->SetRange(3,3);
    TH1 *h_passed = hn_1photon->Projection(1);
    gr[part]->Divide(h_passed, h_total);
    delete h_passed;
    delete h_total;
  }

  mc(0, 2,1);

  gr[0]->SetTitle("Prob Eff for PbSc");
  gr[1]->SetTitle("Prob Eff for PbGl");

  for(int part=0; part<2; part++)
  {
    mcd(0, part+1);
    aset(gr[part], "p_{T} (GeV/c)","Prob Eff", 2.,30., 0.,1.1);
    style(gr[part], part+20, part+1);
    gr[part]->Draw("AP");
  }

  TFile *f_out = new TFile("data/ProbEff-photon.root", "RECREATE");
  for(int part=0; part<2; part++)
    gr[part]->Write();
  f_out->Close();

  c0->Print("plots/ProbEff-photon.pdf");
}
