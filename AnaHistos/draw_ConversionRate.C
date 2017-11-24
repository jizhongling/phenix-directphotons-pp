TGraphErrors* CreateGraph(TFile *f, Int_t criteria)
{
  TH1::SetDefaultSumw2();

  TH2 *h2_conversion = (TH2*)f->Get("h2_conversion");
  TH2 *h2_photon = (TH2*)f->Get("h2_photon");

  const Int_t n = h2_photon->GetNbinsX();
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
    x[i] = h2_photon->GetXaxis()->GetBinCenter(i+1);
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

TGraphErrors* CreateGraph2(TFile *f, Int_t arm)
{
  TH1::SetDefaultSumw2();

  TH2 *h2_conversion = (TH2*)f->Get("h2_conversion");
  TH2 *h2_photon = (TH2*)f->Get("h2_photon");
  TH2 *h2_noconv = (TH2*)f->Get("h2_noconv");

  const Int_t n = h2_photon->GetNbinsX();
  Double_t *x = new Double_t[n];
  Double_t *y = new Double_t[n];
  Double_t *ey = new Double_t[n];
  for(Int_t i=0; i<n; i++)
    x[i] = y[i] = ey[i] = 0.;

  for(Int_t i=0; i<n; i++)
  {
    x[i] = h2_photon->GetXaxis()->GetBinCenter(i+1);
    Double_t conversion1 = h2_conversion->GetBinContent(i+1, arm+1);
    Double_t econversion1 = h2_conversion->GetBinError(i+1, arm+1);
    Double_t conversion2 = h2_conversion->GetBinContent(i+1, arm+3);
    Double_t econversion2 = h2_conversion->GetBinError(i+1, arm+3);
    Double_t photon = h2_photon->GetBinContent(i+1, arm+1);
    Double_t ephoton = h2_photon->GetBinError(i+1, arm+1);
    Double_t noconv = h2_noconv->GetBinContent(i+1, arm+1);
    Double_t enoconv = h2_noconv->GetBinError(i+1, arm+1);
    if( conversion1 > 0. && conversion2 > 0. && photon > 0. && noconv > 0. )
    {
      y[i] = ( conversion1 + conversion2 + noconv ) / photon;
      ey[i] = y[i] * sqrt( pow(econversion1/conversion1,2.) + pow(econversion2/conversion2,2.) + pow(ephoton/photon,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(n, x, y, 0, 0);
  return graph;
}

TGraphErrors* CreateGraph3(TFile *f, Int_t arm)
{
  TH1::SetDefaultSumw2();

  TH2 *h2_photon = (TH2*)f->Get("h2_photon");
  TH2 *h2_vtxconv = (TH2*)f->Get("h2_vtxconv");

  const Int_t n = h2_photon->GetNbinsX();
  Double_t *x = new Double_t[n];
  Double_t *y = new Double_t[n];
  Double_t *ey = new Double_t[n];
  for(Int_t i=0; i<n; i++)
    x[i] = y[i] = ey[i] = 0.;

  for(Int_t i=0; i<n; i++)
  {
    x[i] = h2_photon->GetXaxis()->GetBinCenter(i+1);
    Double_t photon = h2_photon->GetBinContent(i+1, arm+1);
    Double_t ephoton = h2_photon->GetBinError(i+1, arm+1);
    Double_t vtxconv = h2_vtxconv->GetBinContent(i+1, arm+1);
    Double_t evtxconv = h2_vtxconv->GetBinError(i+1, arm+1);
    if( photon > 0. && vtxconv > 0. )
    {
      y[i] = vtxconv / photon;
      ey[i] = y[i] * sqrt( pow(evtxconv/vtxconv,2.) + pow(ephoton/photon,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(n, x, y, 0, ey);
  return graph;
}

void draw_ConversionRate()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/MissingRatio-macros/MissingRatio-histo.root");
  TObjArray *Glist = new TObjArray();

  mc();
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

  legi(0, 0.2, 0.7, 0.8, 0.9);
  //leg0->AddEntry(gr[0], "#splitline{west arm conversion}{large opening angle}", "P");
  //leg0->AddEntry(gr[1], "#splitline{east arm conversion}{large opening angle}", "P");
  //leg0->AddEntry(gr[2], "#splitline{west arm conversion}{small opening angle}", "P");
  //leg0->AddEntry(gr[3], "#splitline{east arm conversion}{small opening angle}", "P");
  leg0->AddEntry(gr[0], "west large opening angle", "P");
  leg0->AddEntry(gr[1], "east large opening angle", "P");
  leg0->AddEntry(gr[2], "west small opening angle", "P");
  leg0->AddEntry(gr[3], "east small opening angle", "P");
  leg0->Draw();

  mc(1);
  mc(2);
  TGraphErrors *gr2[2];
  TGraphErrors *gr3[2];

  for(Int_t arm=0; arm<2; arm++)
  {
    mcd(1);
    gr2[arm] = CreateGraph2(f, arm);
    aset(gr2[arm], "p_{T} [GeV/c]", "rate", 0.,30., 0.8,1.3);
    style(gr2[arm], 20+arm, 1+arm);
    if(arm==0)
      gr2[arm]->Draw("AP");
    else
      gr2[arm]->Draw("P");

    mcd(2);
    gr3[arm] = CreateGraph3(f, arm);
    aset(gr3[arm], "p_{T} [GeV/c]", "rate", 0.,30., 0.,0.3);
    style(gr3[arm], 20+arm, 1+arm);
    if(arm==0)
      gr3[arm]->Draw("AP");
    else
      gr3[arm]->Draw("P");
  }

  mcd(1);
  legi(1, 0.2, 0.7, 0.8, 0.9);
  leg1->AddEntry(gr2[0], "west e^{+}e^{-} pair conv plus no conv", "P");
  leg1->AddEntry(gr2[1], "east e^{+}e^{-} pair conv plus no conv", "P");
  leg1->Draw();

  mcd(2);
  line->DrawLine(0.,0.099,30.,0.099);
  line->DrawLine(0.,0.147,30.,0.147);
  legi(2, 0.2, 0.7, 0.8, 0.9);
  leg2->AddEntry(gr3[0], "west large opening angle", "P");
  leg2->AddEntry(gr3[1], "east large opening angle", "P");
  leg2->Draw();

  c0->Print("ConversionRate.pdf");
  c1->Print("ConversionRate-total.pdf");
  c2->Print("ConversionRate-allconv.pdf");

  TFile *fout = new TFile("data/ConversionRate.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
