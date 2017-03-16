TGraphErrors *DivideHisto(TH1 *h1, TH1 *h2)
{
  Int_t gn = h1->GetXaxis()->GetNbins();
  Double_t *gx = new Double_t[gn];
  Double_t *gy = new Double_t[gn];
  Double_t *egx = new Double_t[gn];
  Double_t *egy = new Double_t[gn];

  for(Int_t i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egx[i] = egy[i] = 0.;

    Double_t h1x = h1->GetXaxis()->GetBinCenter(i+1);
    Double_t h1y = h1->GetBinContent(i+1); 
    Double_t h2y = h2->GetBinContent(i+1); 
    Double_t eh1x = h1->GetXaxis()->GetBinWidth(i+1) / 2.;
    Double_t eh1y = h1->GetBinError(i+1);
    Double_t eh2y = h2->GetBinError(i+1);

    if( h1y > 0. && h2y > 0. )
    {
      gx[i] = h1x;
      gy[i] = h1y / h2y;
      egx[i] = eh1x;
      egy[i] = gy[i] * sqrt( pow(eh1y/h1y,2.) + pow(eh2y/h2y,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(gn, gx, gy, egx, egy);
  return graph;
}

void draw_pTSmear()
{
  TFile *f = new TFile("AnaPHPythia-histo.root");
  TFile *f1 = new TFile("MissingRatio-histo.root");

  TH1::SetDefaultSumw2();

  THnSparse *hn_pi0_total = (THnSparse*)f->Get("hn_total");
  THnSparse *hn_pi0_separate = (THnSparse*)f->Get("hn_separate");
  THnSparse *hn_total = (THnSparse*)f1->Get("hn_total");
  THnSparse *hn_pion = (THnSparse*)f1->Get("hn_pion");

  TCanvas *c = new TCanvas("c", "Canvas", 1200, 600);
  gStyle->SetOptStat(0);
  c->Divide(2,1);

  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};
  const char *pname[2] = {"PbSc", "PbGl"};

  hn_pion->GetAxis(2)->SetRange(112,162);

  for(Int_t part=0; part<2; part++)
  {
    hn_pi0_total->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_pi0_separate->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_total->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_pion->GetAxis(3)->SetRange(secl[part],sech[part]);

    hn_pi0_separate->GetAxis(1)->SetRange(0,31);
    TH1 *h_pi0_truth = hn_pi0_separate->Projection(0)->Clone("h_pi0_truth");
    hn_pi0_separate->GetAxis(0)->SetRange(0,31);
    TH1 *h_pi0_reco = hn_pi0_separate->Projection(1)->Clone("h_pi0_reco");

    hn_pion->GetAxis(1)->SetRange(0,31);
    TH1 *h_truth = hn_pion->Projection(0)->Clone("h_truth");
    hn_pion->GetAxis(0)->SetRange(0,31);
    TH1 *h_reco = hn_pion->Projection(1)->Clone("h_reco");

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

  c->Print("pTSmear.pdf");
}
