void draw_TriggerEfficiency()
{
  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510MinBias/8905/data/total.root");
  TH3 *h3_trig = (TH3*)f->Get("h3_trig");

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TH1 *h_ertbw = h3_trig->ProjectionX("h_ertbw", 1, 4, 5, 5);
  TH1 *h_totalw = h3_trig->ProjectionX("h_totalw", 1, 4, 1, 1);
  TGraphAsymmErrors *gr_ertbw = new TGraphAsymmErrors(h_ertbw, h_totalw);
  gr_ertbw->SetTitle("ERT_4x4b efficiency");
  gr_ertbw->GetXaxis()->SetTitle("p_{T} [GeV]");
  gr_ertbw->GetYaxis()->SetTitle("efficiency");
  gr_ertbw->GetYaxis()->SetTitleOffset(1.2);
  gr_ertbw->GetXaxis()->SetRangeUser(0., 30.);
  gr_ertbw->GetYaxis()->SetRangeUser(0., 1.1);
  gr_ertbw->SetMarkerColor(2);
  gr_ertbw->SetMarkerStyle(20);
  gr_ertbw->SetMarkerSize(1.);
  gr_ertbw->Draw("AP");

  TH1 *h_ertbe = h3_trig->ProjectionX("h_ertbe", 5, 8, 5, 5);
  TH1 *h_totale = h3_trig->ProjectionX("h_totale", 5, 8, 1, 1);
  TGraphAsymmErrors *gr_ertbe = new TGraphAsymmErrors(h_ertbe, h_totale);
  gr_ertbe->SetMarkerColor(4);
  gr_ertbe->SetMarkerStyle(21);
  gr_ertbe->SetMarkerSize(1.);
  gr_ertbe->Draw("P");

  TLegend *leg = new TLegend(0.6, 0.1, 0.9, 0.3);
  leg->AddEntry(gr_ertbw, "West arm", "LPE");
  leg->AddEntry(gr_ertbe, "East arm", "LPE");
  leg->Draw();

  c->Print("TriggerEfficiency.pdf");

  Int_t grnw = gr_ertbw->GetN();
  Double_t *grxw = gr_ertbw->GetX();
  Double_t *gryw = gr_ertbw->GetY();
  cout << "West arm: " << grnw << endl;
  for(Int_t i=0; i<grnw; i++)
    cout << grxw[i] << ",";
  cout << endl;
  for(Int_t i=0; i<grnw; i++)
    cout << gryw[i] << ",";
  cout << endl;

  Int_t grne = gr_ertbe->GetN();
  Double_t *grxe = gr_ertbe->GetX();
  Double_t *grye = gr_ertbe->GetY();
  cout << "East arm: " << grne << endl;
  for(Int_t i=0; i<grne; i++)
    cout << grxe[i] << ",";
  cout << endl;
  for(Int_t i=0; i<grne; i++)
    cout << grye[i] << ",";
  cout << endl;
}
