#include "DivideFunctions.h"

void draw_pTSmear()
{
  TFile *f = new TFile("data/AnaPHPythia-histo.root");
  TFile *f1 = new TFile("data/MissingRatio-histo.root");

  TH1::SetDefaultSumw2();

  THnSparse *hn_pi0_total = (THnSparse*)f->Get("hn_total");
  THnSparse *hn_pi0_separate = (THnSparse*)f->Get("hn_separate");
  THnSparse *hn_total = (THnSparse*)f1->Get("hn_total");
  THnSparse *hn_pion = (THnSparse*)f1->Get("hn_pion");

  TCanvas *c = new TCanvas("c", "Canvas", 1200, 600);
  gStyle->SetOptStat(0);
  c->Divide(2,1);

  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};
  const char *pname[2] = {"PbSc", "PbGl"};

  hn_pion->GetAxis(2)->SetRange(112,162);

  for(int part=0; part<2; part++)
  {
    hn_pi0_total->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_pi0_separate->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_total->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_pion->GetAxis(3)->SetRange(secl[part],sech[part]);

    hn_pi0_separate->GetAxis(1)->SetRange(0,31);
    TH1 *h_pi0_truth = (TH1*)hn_pi0_separate->Projection(0)->Clone("h_pi0_truth");
    hn_pi0_separate->GetAxis(0)->SetRange(0,31);
    TH1 *h_pi0_reco = (TH1*)hn_pi0_separate->Projection(1)->Clone("h_pi0_reco");

    hn_pion->GetAxis(1)->SetRange(0,31);
    TH1 *h_truth = (TH1*)hn_pion->Projection(0)->Clone("h_truth");
    hn_pion->GetAxis(0)->SetRange(0,31);
    TH1 *h_reco = (TH1*)hn_pion->Projection(1)->Clone("h_reco");

    c->cd(part+1);

    TGraphErrors *gr_pi0 = DivideHisto(h_pi0_reco, h_pi0_truth);
    gr_pi0->SetTitle(Form("pT smearing for %s",pname[part]));
    gr_pi0->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    gr_pi0->GetYaxis()->SetTitle("Smearing");
    gr_pi0->GetXaxis()->SetRangeUser(0., 30.);
    gr_pi0->GetYaxis()->SetRangeUser(0.8, 1.5);
    gr_pi0->SetMarkerStyle(20);
    gr_pi0->SetMarkerColor(1);
    gr_pi0->Draw("AP");

    TGraphErrors *gr_pisa = DivideHisto(h_reco, h_truth);
    gr_pisa->SetMarkerStyle(21);
    gr_pisa->SetMarkerColor(2);
    gr_pisa->Draw("P");

    TLegend *leg = new TLegend(0.1,0.7,0.4,0.9);
    leg->AddEntry(gr_pi0, "FastMC", "P");
    leg->AddEntry(gr_pisa, "PISA", "P");
    leg->Draw();
  }

  c->Print("plots/pTSmear.pdf");
}
