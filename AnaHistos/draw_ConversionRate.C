TGraphErrors* CreateGraph(TFile *f, Int_t criteria)
{
  TH1::SetDefaultSumw2();

  TH2 *h2_conversion = (TH2*)f->Get("h2_conversion");
  TH2 *h2_photon = (TH2*)f->Get("h2_photon");

  const Int_t n = h2_conversion->GetNbinsX();
  Double_t *x = new Double_t[n];
  Double_t *y = new Double_t[n];
  Double_t *ey = new Double_t[n];
  for(Int_t i=0; i<n; i++)
    x[i] = y[i] = ey[i] = 0.;

  cout << "\nConvertion rate for criteria " << criteria << ": ";
  for(Int_t i=0; i<n; i++)
  {
    x[i] = h2_conversion->GetXaxis()->GetBinCenter(i+1);
    Double_t conversion = h2_conversion->GetBinContent(i+1, criteria+1);
    Double_t econversion = h2_conversion->GetBinError(i+1, criteria+1);
    Double_t photon = h2_photon->GetBinContent(i+1, criteria+1);
    Double_t ephoton = h2_photon->GetBinError(i+1, criteria+1);
    if( conversion > 0. && photon > 0. )
    {
      y[i] = conversion / photon;
      ey[i] = y[i] * sqrt( pow(econversion/conversion,2.) + pow(ephoton/photon,2.) );
    }
    cout << y[i] << ", ";
  }
  cout << endl;

  TGraphErrors *graph = new TGraphErrors(n, x, y, 0, ey);
  return graph;
}

void draw_ConversionRate()
{
  TFile *f = new TFile("MissingRatio-histo.root");

  mc();
  legi(0, 0.2, 0.5, 0.8, 0.9);

  TGraphErrors *gr[4];
  for(Int_t icr=0; icr<4; icr++)
  {
    gr[icr] = CreateGraph(f, icr);
    mcd();
    aset(gr[icr], "p_{T} [GeV/c]", "rate", 0.,30., 0.,0.1);
    gr[icr]->SetTitle("Photon conversion rate");
    style(gr[icr], 20+icr, 1+icr);
    if(icr==0)
      gr[icr]->Draw("AP");
    else
      gr[icr]->Draw("P");
  }

  leg0->AddEntry(gr[0], "#splitline{west arm conversion}{outside magnetic field}", "P");
  leg0->AddEntry(gr[1], "#splitline{west arm total conversion}{}", "P");
  leg0->AddEntry(gr[2], "#splitline{east arm conversion}{outside magnetic field}", "P");
  leg0->AddEntry(gr[3], "#splitline{east arm total conversion}{}", "P");
  leg0->Draw();

  c0->Print("ConversionRate.pdf");
}
