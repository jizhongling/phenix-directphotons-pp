#include "GlobalVars.h"

void draw_ToFEff_Photon()
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};
  const char *name[2] = {"PbSc", "PbGl"};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  TAxis *axis_sec = hn_1photon->GetAxis(0);
  TAxis *axis_pt = hn_1photon->GetAxis(1);
  TAxis *axis_pattern = hn_1photon->GetAxis(2);
  TAxis *axis_cut = hn_1photon->GetAxis(3);
  TAxis *axis_type = hn_1photon->GetAxis(4);

  TGraphAsymmErrors *gr[2];
  for(Int_t part=0; part<2; part++)
    gr[part] =  new TGraphAsymmErrors(npT);

  for(Int_t part=0; part<2; part++)
  {
    axis_type->SetRange(3,3);
    axis_sec->SetRange(secl[part],sech[part]);
    axis_cut->SetRange(3,3);
    TH1 *h_total = hn_1photon->Projection(1);
    axis_cut->SetRange(4,4);
    TH1 *h_passed = hn_1photon->Projection(1);
    gr[part]->Divide(h_passed, h_total);
    delete h_passed;
    delete h_total;
  }

  mc(0, 2,1);

  for(Int_t part=0; part<2; part++)
  {
    mcd(0, part+1);
    gr[part]->SetTitle( Form("ToF efficeincy for %s", name[part]) );
    aset(gr[part], "p_{T} [GeV]", "Eff", 2.,30., 0.,1.1);
    style(gr[part], 24, kRed);
    gr[part]->Draw("APE");
    //gr[part]->Fit("pol0", "Q","", 9.,30.);

    //gPad->Update();
    //TPaveStats *st = (TPaveStats*)gr[part]->FindObject("stats");
    //st->SetY1NDC(0.4);
    //st->SetY2NDC(0.6);
  }

  c0->Print("plots/ToFEff-photon.pdf");
}
