void draw_Acceptance()
{
  TH1 *h_phsig[3];
  TFile *fsig = new TFile("Acceptance-signal.root");
  THnSparse *hn_phsig = (THnSparse*)fsig->Get("hn_photon");
  hn_phsig->GetAxis(1)->SetRange(1,4);
  h_phsig[0] = (TH1*)hn_phsig->Projection(0)->Clone("h_phsig_PbScW");
  hn_phsig->GetAxis(1)->SetRange(5,6);
  h_phsig[1] = (TH1*)hn_phsig->Projection(0)->Clone("h_phsig_PbScE");
  hn_phsig->GetAxis(1)->SetRange(7,8);
  h_phsig[2] = (TH1*)hn_phsig->Projection(0)->Clone("h_phsig_PbGlE");
  for(Int_t i=0; i<3; i++)
    h_phsig[i]->SetBinContent(31, 0.);

  TFile *ftot = new TFile("Acceptance-total.root");
  TH1 *h_phtot = (TH1*)ftot->Get("h_photon");
  h_phtot->SetBinContent(31, 1.);

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TGraphAsymmErrors *gr[3];
  for(Int_t i=0; i<3; i++)
  {
    gr[i] = new TGraphAsymmErrors(h_phsig[i], h_phtot);
    gr[i]->SetTitle("Acceptance");
    gr[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr[i]->GetYaxis()->SetTitle("Acceptance");
    gr[i]->GetYaxis()->SetTitleOffset(1.2);
    gr[i]->GetXaxis()->SetRangeUser(0., 30.);
    gr[i]->GetYaxis()->SetRangeUser(0., 0.3);
    gr[i]->SetMarkerColor(i+1);
    gr[i]->SetMarkerStyle(i+20);
    if(i==0)
      gr[i]->Draw("AP");
    else
      gr[i]->Draw("P");
  }

  TLegend *leg = new TLegend(0.4, 0.7, 0.6, 0.9);
  leg->AddEntry(gr[0], "PbScW", "LEP");
  leg->AddEntry(gr[1], "PbScE", "LEP");
  leg->AddEntry(gr[2], "PbGlE", "LEP");
  leg->Draw();

  c->Print("Acceptance.pdf");

  Int_t n[3];
  Double_t *x[3];
  Double_t *y[3];
  Double_t *eyhigh[3];
  Double_t *eylow[3];
  for(Int_t i=0; i<3; i++)
  {
    n[i] = gr[i]->GetN();
    x[i] = gr[i]->GetX();
    y[i] = gr[i]->GetY();
    eyhigh[i] = gr[i]->GetEYhigh();
    eylow[i] = gr[i]->GetEYlow();
  }

  for(Int_t i=0; i<3; i++)
  {
    cout << "\nPart " << i << " Y: ";
    for(Int_t j=0; j<n[i]; j++)
      cout << y[i][j] << ",";
  }
  cout << endl;

  for(Int_t i=0; i<3; i++)
  {
    cout << "\nPart " << i << " EY: ";
    for(Int_t j=0; j<n[i]; j++)
    {
      Double_t ey = eyhigh[i][j] > eylow[i][j] ? eyhigh[i][j] : eylow[i][j];
      cout << ey << ",";
    }
  }
  cout << endl;
}
