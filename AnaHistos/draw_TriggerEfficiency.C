void GenerateTriggerEfficiency(TFile *f, Int_t ispion)
{
  TH3 *h3_trig;
  if(ispion == 0)
    h3_trig = (TH3*)f->Get("h3_trig");
  else if(ispion == 1)
    h3_trig = (TH3*)f->Get("h3_trig_pion");

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TH1 *h_ertb_PbScW = h3_trig->ProjectionX("h_ertb_PbScW", 1, 4, 9, 9);
  TH1 *h_total_PbScW = h3_trig->ProjectionX("h_total_PbScW", 1, 4, 1, 1);
  TGraphAsymmErrors *gr_ertb_PbScW = new TGraphAsymmErrors(h_ertb_PbScW, h_total_PbScW);
  gr_ertb_PbScW->SetTitle("ERT_4x4b efficiency");
  if(ispion == 0)
    gr_ertb_PbScW->SetTitle("Trigger efficiency for photon");
  else if(ispion == 1)
    gr_ertb_PbScW->SetTitle("Trigger efficeincy for #pi^{0}");
  gr_ertb_PbScW->GetXaxis()->SetTitle("p_{T} [GeV]");
  gr_ertb_PbScW->GetYaxis()->SetTitle("efficiency");
  gr_ertb_PbScW->GetYaxis()->SetTitleOffset(1.2);
  gr_ertb_PbScW->GetXaxis()->SetRangeUser(0., 30.);
  gr_ertb_PbScW->GetYaxis()->SetRangeUser(0., 1.1);
  gr_ertb_PbScW->SetMarkerColor(1);
  gr_ertb_PbScW->SetMarkerStyle(20);
  gr_ertb_PbScW->SetMarkerSize(1.);
  gr_ertb_PbScW->Draw("AP");

  TH1 *h_ertb_PbScE = h3_trig->ProjectionX("h_ertb_PbScE", 5, 6, 9, 9);
  TH1 *h_total_PbScE = h3_trig->ProjectionX("h_total_PbScE", 5, 6, 1, 1);
  TGraphAsymmErrors *gr_ertb_PbScE = new TGraphAsymmErrors(h_ertb_PbScE, h_total_PbScE);
  gr_ertb_PbScE->SetMarkerColor(2);
  gr_ertb_PbScE->SetMarkerStyle(21);
  gr_ertb_PbScE->SetMarkerSize(1.);
  gr_ertb_PbScE->Draw("P");

  TH1 *h_ertb_PbGlE = h3_trig->ProjectionX("h_ertb_PbGlE", 7, 8, 9, 9);
  TH1 *h_total_PbGlE = h3_trig->ProjectionX("h_total_PbGlE", 7, 8, 1, 1);
  TGraphAsymmErrors *gr_ertb_PbGlE = new TGraphAsymmErrors(h_ertb_PbGlE, h_total_PbGlE);
  gr_ertb_PbGlE->SetMarkerColor(3);
  gr_ertb_PbGlE->SetMarkerStyle(22);
  gr_ertb_PbGlE->SetMarkerSize(1.);
  gr_ertb_PbGlE->Draw("P");

  TLegend *leg = new TLegend(0.6, 0.1, 0.9, 0.3);
  leg->AddEntry(gr_ertb_PbScW, "PbScW", "LPE");
  leg->AddEntry(gr_ertb_PbScE, "PbScE", "LPE");
  leg->AddEntry(gr_ertb_PbGlE, "PbGlE", "LPE");
  leg->Draw();

  if(ispion == 0)
  {
    c->Print("TriggerEfficiency-photon.pdf");
    delete c;
  }
  else if(ispion == 1)
  {
    c->Print("TriggerEfficiency-pion.pdf");
    delete c;
  }

  Int_t grn_PbScW = gr_ertb_PbScW->GetN();
  Double_t *grx_PbScW = gr_ertb_PbScW->GetX();
  Double_t *gry_PbScW = gr_ertb_PbScW->GetY();
  Double_t *egry_PbScWhigh = gr_ertb_PbScW->GetEYhigh();
  Double_t *egry_PbScWlow = gr_ertb_PbScW->GetEYlow();
  cout << "PbScW: " << grn_PbScW;
  cout << "\nX: ";
  for(Int_t i=0; i<grn_PbScW; i++)
    cout << grx_PbScW[i] << ",";
  cout << "\nY: ";
  for(Int_t i=0; i<grn_PbScW; i++)
    cout << gry_PbScW[i] << ",";
  cout << "\nEY: ";
  for(Int_t i=0; i<grn_PbScW; i++)
    cout << ( egry_PbScWhigh[i] > egry_PbScWlow[i] ? egry_PbScWhigh[i] : egry_PbScWlow[i] ) << ",";
  cout << endl;

  Int_t grn_PbScE = gr_ertb_PbScE->GetN();
  Double_t *grx_PbScE = gr_ertb_PbScE->GetX();
  Double_t *gry_PbScE = gr_ertb_PbScE->GetY();
  Double_t *egry_PbScEhigh = gr_ertb_PbScE->GetEYhigh();
  Double_t *egry_PbScElow = gr_ertb_PbScE->GetEYlow();

  cout << "\nPbScE: " << grn_PbScE;
  cout << "\nX: ";
  for(Int_t i=0; i<grn_PbScE; i++)
    cout << grx_PbScE[i] << ",";
  cout << "\nY: ";
  for(Int_t i=0; i<grn_PbScE; i++)
    cout << gry_PbScE[i] << ",";
  cout << "\nEY: ";
  for(Int_t i=0; i<grn_PbScE; i++)
    cout << ( egry_PbScEhigh[i] > egry_PbScElow[i] ? egry_PbScEhigh[i] : egry_PbScElow[i] ) << ",";
  cout << endl;

  Int_t grn_PbGlE = gr_ertb_PbGlE->GetN();
  Double_t *grx_PbGlE = gr_ertb_PbGlE->GetX();
  Double_t *gry_PbGlE = gr_ertb_PbGlE->GetY();
  Double_t *egry_PbGlEhigh = gr_ertb_PbGlE->GetEYhigh();
  Double_t *egry_PbGlElow = gr_ertb_PbGlE->GetEYlow();
  cout << "\nPbGlE: " << grn_PbGlE;
  cout << "\nX: ";
  for(Int_t i=0; i<grn_PbGlE; i++)
    cout << grx_PbGlE[i] << ",";
  cout << "\nY: ";
  for(Int_t i=0; i<grn_PbGlE; i++)
    cout << gry_PbGlE[i] << ",";
  cout << "\nEY: ";
  for(Int_t i=0; i<grn_PbGlE; i++)
    cout << ( egry_PbGlEhigh[i] > egry_PbGlElow[i] ? egry_PbGlEhigh[i] : egry_PbGlElow[i] ) << ",";
  cout << endl;

  return;
}

void draw_TriggerEfficiency()
{
  //TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510MinBias/9161/data/total.root");
  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/9473/data/total.root");
  for(Int_t ispion=0; ispion<2; ispion++)
  {
    cout << "\nispion " << ispion << endl;
    GenerateTriggerEfficiency(f, ispion);
  }
}
