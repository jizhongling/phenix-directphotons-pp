#include "GlobalVars.h"

void draw_ERTEff_Photon()
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

  TH3 *h3_ert = (TH3*)f->Get("h3_ert");

  for(int part=0; part<2; part++)
  {
    TH1 *h_total = (TH1*)h3_ert->ProjectionY("_py", secl[part],sech[part], 3,3)->Clone("h_total");
    TH1 *h_passed = (TH1*)h3_ert->ProjectionY("_py", secl[part],sech[part], 6,6)->Clone("h_passed");
    gr[part]->Divide(h_passed, h_total);
  }

  mc();
  mcd();

  for(int part=0; part<2; part++)
  {
    gr[part]->SetName(Form("gr_%d",part));
    gr[part]->SetTitle("ERT_4x4c trigger efficeincy for photon");
    aset(gr[part], "p_{T} [GeV]","Eff", 1.,30., 0.,1.1);
    style(gr[part], part+20, part+1);
    if(part==0)
      gr[part]->Draw("APE");
    else
      gr[part]->Draw("PE");

    gr[part]->Fit("pol0", "Q","", 10.,30.);
    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr[part]->FindObject("stats");
    st->SetY1NDC(0.3-part*0.2);
    st->SetY2NDC(0.5-part*0.2);
  }

  legi(0, 0.7,0.2,0.9,0.4);
  leg0->AddEntry(gr[0], "PbSc", "LPE");
  leg0->AddEntry(gr[1], "PbGl", "LPE");
  leg0->Draw();

  c0->Print("plots/ERTEff-photon.pdf");

  TFile *f_out = new TFile("data/ERTEff-photon.root", "RECREATE");
  for(int part=0; part<2; part++)
    gr[part]->Write();
  f_out->Close();
}
