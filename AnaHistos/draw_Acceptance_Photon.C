#include "GlobalVars.h"

void draw_Acceptance_Photon()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TGraphAsymmErrors *gr[3];
  for(int part=0; part<3; part++)
  {
    gr[part] = new TGraphAsymmErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");

  TH1 *h_total = (TH1*)f->Get("h_photon");
  THnSparse *hn_photon = (THnSparse*)f->Get("hn_photon");

  for(int part=0; part<3; part++)
  {
    hn_photon->GetAxis(2)->SetRange(secl[part],sech[part]);
    TH1 *h_pass = hn_photon->Projection(1);
    gr[part]->Divide(h_pass, h_total, "n");
    delete h_pass;
  }

  mc();
  mcd();
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(int part=0; part<3; part++)
  {
    gr[part]->SetTitle("Photon acceptance");
    aset(gr[part], "p_{T} [GeV]","acceptance", 2.,30., 0.,0.2);
    style(gr[part], part+20, part+1);
    if(part==0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
    leg0->AddEntry(gr[part], pname[part], "P");
  }

  leg0->Draw();
  c0->Print("plots/Acceptance-photon.pdf");

  TFile *f_out = new TFile("data/Acceptance-photon.root", "RECREATE");
  for(int part=0; part<3; part++)
  {
    gr[part]->Write();
  }
  f_out->Close();
}
