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

  Double_t tot_conv = 0.;
  Double_t tot_ph = 0.;
  cout << "\nConvertion rate for criteria " << criteria << ": ";
  for(Int_t i=0; i<n; i++)
  {
    x[i] = h2_conversion->GetXaxis()->GetBinCenter(i+1);
    Double_t conversion = h2_conversion->GetBinContent(i+1, criteria+1);
    Double_t econversion = h2_conversion->GetBinError(i+1, criteria+1);
    Double_t photon = h2_photon->GetBinContent(i+1, criteria+1);
    Double_t ephoton = h2_photon->GetBinError(i+1, criteria+1);
    if(i>25)
    {
      tot_conv += conversion;
      tot_ph += photon;
    }
    if( conversion > 0. && photon > 0. )
    {
      y[i] = conversion / photon;
      ey[i] = y[i] * sqrt( pow(econversion/conversion,2.) + pow(ephoton/photon,2.) );
    }
    cout << y[i] << ", ";
  }
  cout << endl << "Criteria = " << criteria << "\tRate = " << tot_conv/tot_ph << endl;

  TGraphErrors *graph = new TGraphErrors(n, x, y, 0, ey);
  return graph;
}

void draw_ConversionRate()
{
  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");
  TObjArray *Glist = new TObjArray();

  mc();
  legi(0, 0.2, 0.7, 0.8, 0.9);

  TGraphErrors *gr[4];
  for(Int_t icr=0; icr<4; icr++)
  {
    gr[icr] = CreateGraph(f, icr);
    Glist->Add(gr[icr]);
    gr[icr]->SetName(Form("gr_%d",icr));
    mcd();
    aset(gr[icr], "p_{T} [GeV/c]", "rate", 0.,30., 0.,0.3);
    gr[icr]->SetTitle("Photon conversion rate");
    style(gr[icr], 20+icr, 1+icr);
    if(icr==0)
      gr[icr]->Draw("AP");
    else
      gr[icr]->Draw("P");
  }

  TLine *line = new TLine();
  line->SetLineColor(kRed);
  line->SetLineWidth(5);
  line->DrawLine(0.,0.099,30.,0.099);
  line->DrawLine(0.,0.147,30.,0.147);

  //leg0->AddEntry(gr[0], "#splitline{west arm conversion}{large opening angle}", "P");
  //leg0->AddEntry(gr[1], "#splitline{east arm conversion}{large opening angle}", "P");
  //leg0->AddEntry(gr[2], "#splitline{west arm conversion}{small opening angle}", "P");
  //leg0->AddEntry(gr[3], "#splitline{east arm conversion}{small opening angle}", "P");
  leg0->AddEntry(gr[0], "west large opening angle", "P");
  leg0->AddEntry(gr[1], "east large opening angle", "P");
  leg0->AddEntry(gr[2], "west small opening angle", "P");
  leg0->AddEntry(gr[3], "east small opening angle", "P");
  leg0->Draw();

  c0->Print("ConversionRate.pdf");

  TFile *fout = new TFile("ConversionRate.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
