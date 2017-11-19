#include "DivideFunctions.h"

void draw_Merge()
{
  const Int_t npT = 31;
  const Double_t vpT[npT] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

  cross = new TF1("cross", "x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))", 0, 30);
  cross->SetParameters(2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00);

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");
  TFile *f1 = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");
  //TFile *f1 = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");

  const Int_t reco = 0; // truth pT:0; reconstructed pT: 1

  TH1::SetDefaultSumw2();

  THnSparse *hn_pi0_total = (THnSparse*)f->Get("hn_total");
  THnSparse *hn_pi0_separate = (THnSparse*)f->Get("hn_separate");
  THnSparse *hn_pi0_pion = (THnSparse*)f->Get("hn_pion");

  THnSparse *hn_total = (THnSparse*)f1->Get("hn_total");
  THnSparse *hn_separate = (THnSparse*)f1->Get("hn_separate");
  THnSparse *hn_pion = (THnSparse*)f1->Get("hn_pion");

  mc();
  legi(0, 0.2, 0.2, 0.5, 0.5);

  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};
  const char *pname[2] = {"PbSc", "PbGl"};

  for(Int_t part=0; part<2; part++)
  {
    hn_pi0_total->GetAxis(2)->SetRange(secl[part], sech[part]);
    TH1 *h_pi0_total = hn_pi0_total->Projection(reco);
    hn_pi0_separate->GetAxis(2)->SetRange(secl[part], sech[part]);
    TH1 *h_pi0_separate = hn_pi0_separate->Projection(reco);
    //hn_pi0_pion->GetAxis(3)->SetRange(secl[part], sech[part]);
    //TH1 *h_pi0_separate = new TH1F("h_pi0_separate", "Reconstructed #pi^{0}", npT-1, vpT);
    //for(Int_t ipt=0; ipt<30; ipt++)
    //{
    //  hn_pi0_pion->GetAxis(reco)->SetRange(ipt+1,ipt+1);
    //  TH1 *h_minv = hn_pi0_pion->Projection(2);
    //  Double_t e1, e2, e3;
    //  Double_t n_separate = h_minv->IntegralAndError(113,162,e1) - ( h_minv->IntegralAndError(48,97,e2) - h_minv->IntegralAndError(178,227,e3) ) / 2.;
    //  Double_t e_separate = sqrt( e1*e1 + e2*e2/4. + e3*e3/4. );
    //  h_pi0_separate->SetBinContent(ipt+1, n_separate);
    //  h_pi0_separate->SetBinError(ipt+1, e_separate);
    //  delete h_minv;
    //}

    //TGraphAsymmErrors *gr_pi0 = new TGraphAsymmErrors(h_pi0_separate, h_pi0_total, "n");
    TGraphErrors *gr_pi0 = DivideHisto(h_pi0_separate, h_pi0_total);
    //InvertGraph(gr_pi0);
    delete h_pi0_separate;
    delete h_pi0_total;

    hn_total->GetAxis(2)->SetRange(secl[part], sech[part]);
    TH1 *h_total = hn_total->Projection(reco);
    hn_separate->GetAxis(2)->SetRange(secl[part], sech[part]);
    TH1 *h_separate = hn_separate->Projection(reco);
    //hn_pion->GetAxis(3)->SetRange(secl[part], sech[part]);
    //TH1 *h_separate = new TH1F("h_separate", "Reconstructed #pi^{0}", npT-1, vpT);
    //for(Int_t ipt=0; ipt<30; ipt++)
    //{
    //  hn_pion->GetAxis(reco)->SetRange(ipt+1,ipt+1);
    //  TH1 *h_minv = hn_pion->Projection(2);
    //  Double_t e1, e2, e3;
    //  Double_t n_separate = h_minv->IntegralAndError(0,301,e1); // - ( h_minv->IntegralAndError(48,97,e2) - h_minv->IntegralAndError(178,227,e3) ) / 2.;
    //  Double_t e_separate = sqrt( e1*e1 ); // + e2*e2/4. + e3*e3/4. );
    //  h_separate->SetBinContent(ipt+1, n_separate);
    //  h_separate->SetBinError(ipt+1, e_separate);
    //  delete h_minv;
    //}

    //TGraphAsymmErrors *gr_pion_1 = new TGraphAsymmErrors(h_separate, h_total, "n");
    TGraphErrors *gr_pion_1 = DivideHisto(h_separate, h_total);
    //InvertGraph(gr_pion_1);
    delete h_separate;
    delete h_total;

    mcd();
    aset(gr_pi0, "p_{T} [GeV/c]", "rate");
    style(gr_pi0, 20+part, 1+part);
    gr_pi0->SetTitle("Separating rate for #pi^{0}");
    gr_pi0->GetXaxis()->SetRangeUser(0., 30.);
    gr_pi0->GetYaxis()->SetRangeUser(0., 1.);
    if(part==0)
      gr_pi0->Draw("AP");
    else
      gr_pi0->Draw("P");

    style(gr_pion_1, 24+part, 1+part);
    //gr_pion_1->Draw("P");

    leg0->AddEntry(gr_pi0, Form("%s", pname[part]), "P");
    //leg0->AddEntry(gr_pion_1, Form("PISA, %s", pname[part]), "P");
    leg0->Draw();
  }

  c0->Print("Merge.pdf");
}
