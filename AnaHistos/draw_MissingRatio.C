void GetCriteria(Int_t criteria, Int_t &merge, Int_t &part, Int_t &prob, Int_t &warnmap)
{
  merge = criteria / 12;
  part = criteria % 12 / 4;
  prob = criteria % 4 / 2;
  warnmap = criteria % 2;

  return;
}

TGraphErrors* CreateGraph(TFile *f, const Int_t criteria, Int_t ispion)
{
  TH1::SetDefaultSumw2();

  if(ispion == 0)
  {
    h2_missing = (TH2*)f->Get("h2_missing");
    h2_nomissing = (TH2*)f->Get("h2_nomissing");
  }
  else if(ispion == 1)
  {
    h2_missing = (TH2*)f->Get("h2_missing_pion");
    h2_nomissing = (TH2*)f->Get("h2_nomissing_pion");
  }

  const Int_t n = h2_missing->GetNbinsX();
  Double_t *x = new Double_t[n];
  Double_t *y = new Double_t[n];
  Double_t *ey = new Double_t[n];
  for(Int_t i=0; i<n; i++)
    x[i] = y[i] = ey[i] = 0.;

  Int_t merge, part, prob, warnmap;
  GetCriteria(criteria, merge, part, prob, warnmap);

  for(Int_t i=0; i<n; i++)
  {
    x[i] = h2_missing->GetXaxis()->GetBinCenter(i+1);
    Double_t missing = h2_missing->GetBinContent(i+1, criteria+1);
    Double_t emissing = h2_missing->GetBinError(i+1, criteria+1);
    Double_t nomissing = h2_nomissing->GetBinContent(i+1, criteria+1);
    Double_t enomissing = h2_nomissing->GetBinError(i+1, criteria+1);
    if( missing > 0. && nomissing > 0. )
    {
      y[i] = missing / nomissing;
      ey[i] = y[i] * sqrt( pow(emissing/missing,2.) + pow(enomissing/nomissing,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(n, x, y, 0, ey);
  return graph;
}

void GenerateMissingRatio(TFile *f, TObjArray *Glist, Int_t ispion)
{
  TCanvas *c = new TCanvas("c", "Canvas", 1800, 1200);
  gStyle->SetOptStat(0);
  c->Divide(3,2);

  TGraphErrors *gr[24];
  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);

  for(Int_t icr=0; icr<24; icr++)
  {
    c->cd( icr/4 + 1 );

    char buf[100];
    Int_t merge, part, prob, warnmap;
    GetCriteria(icr, merge, part, prob, warnmap);

    gr[icr] = CreateGraph(f, icr, ispion);
    Glist->Add(gr[icr]);
    gr[icr]->SetName(Form("gr_%d",24*ispion+icr));
    sprintf(buf, "Missing ratio for merge %d part %d", merge, part);
    gr[icr]->SetTitle(buf);
    gr[icr]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    gr[icr]->GetYaxis()->SetTitle("ratio");
    gr[icr]->GetXaxis()->SetRangeUser(0., 30.);
    gr[icr]->GetYaxis()->SetRangeUser(0., 2.);
    gr[icr]->SetMarkerColor( icr%4 + 1 );
    gr[icr]->SetMarkerStyle( icr%4 + 20 );
    gr[icr]->SetMarkerSize(1.5);
    if( icr%4 == 0 )
      gr[icr]->Draw("AP");
    else
      gr[icr]->Draw("P");

    sprintf(buf, "prob %d warnmap %d", prob, warnmap);
    leg->AddEntry(gr[icr], buf, "P");
    if( icr%4 == 3 )
    {
      leg->DrawClone();
      leg->Clear();
    }
  }

  if(ispion == 0)
  {
    c->Print("plots/MissingRatio-photon.pdf");
    delete c;
  }
  else if(ispion == 1)
  {
    c->Print("plots/MissingRatio-pion.pdf");
    delete c;
  }

  return;
}

void draw_MissingRatio()
{
  TFile *f = new TFile("data/MissingRatio-histo.root");
  //TFile *f = new TFile("data/Acceptance-histo.root");
  //TFile *f = new TFile("data/AnaPHPythia-histo.root");
  TObjArray *Glist = new TObjArray();

  for(Int_t ispion=0; ispion<2; ispion++)
  {
    cout << "\nispion " << ispion << endl;
    GenerateMissingRatio(f, Glist, ispion);
  }

  TFile *fout = new TFile("data/MissingRatio.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
