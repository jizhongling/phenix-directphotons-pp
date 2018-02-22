#include "GlobalVars.h"

void draw_BBCEff_Photon()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TGraphAsymmErrors *gr[2];
  for(Int_t part=0; part<2; part++)
  {
    gr[part] = new TGraphAsymmErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  TH3 *h3_trig = (TH3*)f->Get("h3_bbc");

  for(Int_t part=0; part<2; part++)
  {
    TH1 *h_total = (TH1*)h3_trig->ProjectionY("_py", secl[part],sech[part], 1,1)->Clone("h_total");
    TH1 *h_passed = (TH1*)h3_trig->ProjectionY("_py", secl[part],sech[part], 2,2)->Clone("h_passed");
    gr[part]->Divide(h_passed, h_total);
    delete h_total;
    delete h_passed;
  }

  mc(0, 2,1);

  for(Int_t part=0; part<2; part++)
  {
    mcd(0, part+1);
    gr[part]->SetTitle( Form("BBC trigger efficeincy for %s",pname[part]) );
    aset(gr[part], "p_{T} [GeV]","Eff", 0.,30., 0.,1.);
    style(gr[part], part+20, part+1);
    gr[part]->Draw("APE");
    //gr[part]->Fit("pol0", "Q","", 2.,20.);

    //gPad->Update();
    //TPaveStats *st = (TPaveStats*)gr[part]->FindObject("stats");
    //st->SetY1NDC(0.6);
    //st->SetY2NDC(0.8);
  }

  TFile *f_out = new TFile("data/BBCEff-photon.root", "RECREATE");
  for(Int_t part=0; part<2; part++)
    gr[part]->Write();
  f_out->Close();

  c0->Print("plots/BBCEff-photon.pdf");
}
