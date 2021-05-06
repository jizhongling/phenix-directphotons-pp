#include "DivideFunctions.h"

void draw_ProcessRatio()
{
  gSystem->Load("libTOAD.so");

  const int proc_qcd[6] = {11, 12, 13, 28, 53, 68};
  const int proc_photon[5] = {14, 18, 29, 114, 115};

  TOAD *toad_loader = new TOAD("AnaFastMC");
  string file_weights = toad_loader->location("MinBiasPtWeights.root");

  TFile *f = new TFile( file_weights.c_str() );
  TH2 *h2_proc_pt = (TH2*)f->Get("h2_proc_pt");

  TH1 *h_total = (TH1*)h2_proc_pt->ProjectionX()->Clone("h_total");
  h_total->Rebin(10);
  TH1 *h_qcd = (TH1*)h_total->Clone("h_qcd");
  TH1 *h_photon = (TH1*)h_total->Clone("h_photon");
  h_qcd->Reset();
  h_photon->Reset();

  for(int i=0; i<6; i++)
  {
    TH1 *h_tmp = h2_proc_pt->ProjectionX("h_tmp", proc_qcd[i],proc_qcd[i]);
    h_tmp->Rebin(10);
    h_qcd->Add(h_tmp);
    delete h_tmp;
  }

  for(int i=0; i<5; i++)
  {
    TH1 *h_tmp = h2_proc_pt->ProjectionX("h_tmp", proc_photon[i],proc_photon[i]);
    h_tmp->Rebin(10);
    h_photon->Add(h_tmp);
    delete h_tmp;
  }

  mc();
  mcd();
  TGraphErrors *gr = DivideHisto(h_photon, h_total);
  gr->SetTitle("#frac{Prompt photons}{Hard QCD}");
  aset(gr, "p_{T} (GeV/c)","Ratio", 5.1,30., 0.,0.002);
  style(gr, 20, 1);
  gr->Draw("APE");
  c0->Print("plots/ProcessRatio.pdf");
}
