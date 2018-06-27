void GetCriteria(int criteria, int &merge, int &part, int &prob, int &warnmap)
{
  merge = criteria / 12;
  part = criteria % 12 / 4;
  prob = criteria % 4 / 2;
  warnmap = criteria % 2;

  return;
}

TGraphErrors* CreateGraph(TFile *f, const int criteria, int ispion)
{
  TH1::SetDefaultSumw2();

  if(ispion==0)
  {
    TH2 *h2_missing = (TH2*)f->Get("h2_measured");
    TH2 *h2_nomissing = (TH2*)f->Get("h2_incident");
  }
  else if(ispion==1)
  {
    TH2 *h2_missing = (TH2*)f->Get("h2_measured_pion");
    TH2 *h2_nomissing = (TH2*)f->Get("h2_incident_pion");
  }

  const int n = h2_missing->GetNbinsX();
  double *x = new double[n];
  double *y = new double[n];
  double *ey = new double[n];
  for(int i=0; i<n; i++)
    x[i] = y[i] = ey[i] = 0.;

  int merge, part, prob, warnmap;
  GetCriteria(criteria, merge, part, prob, warnmap);

  for(int i=0; i<n; i++)
  {
    x[i] = h2_missing->GetXaxis()->GetBinCenter(i+1);
    double missing = h2_missing->GetBinContent(i+1, criteria+1);
    double emissing = h2_missing->GetBinError(i+1, criteria+1);
    double nomissing = h2_nomissing->GetBinContent(i+1, criteria+1);
    double enomissing = h2_nomissing->GetBinError(i+1, criteria+1);
    if( missing > 0. && nomissing > 0. )
    {
      y[i] = missing / nomissing;
      ey[i] = y[i] * sqrt( pow(emissing/missing,2.) + pow(enomissing/nomissing,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(n, x, y, 0, ey);
  return graph;
}

void GenerateSmear(TFile *f, TObjArray *Glist, int ispion)
{
  TCanvas *c = new TCanvas("c", "Canvas", 1800, 1200);
  gStyle->SetOptStat(0);
  c->Divide(3,2);

  TGraphErrors *gr[24];
  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);

  for(int icr=0; icr<24; icr++)
  {
    c->cd( icr/4 + 1 );

    char buf[100];
    int merge, part, prob, warnmap;
    GetCriteria(icr, merge, part, prob, warnmap);

    gr[icr] = CreateGraph(f, icr, ispion);
    Glist->Add(gr[icr]);
    gr[icr]->SetName(Form("gr_%d",24*ispion+icr));
    sprintf(buf, "Smearing factor for merge %d part %d", merge, part);
    gr[icr]->SetTitle(buf);
    gr[icr]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    gr[icr]->GetYaxis()->SetTitle("#epsilon_{smear}");
    gr[icr]->GetXaxis()->SetRangeUser(0., 30.);
    gr[icr]->GetYaxis()->SetRangeUser(0.5, 1.5);
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

  if(ispion==0)
  {
    c->Print("plots/Smear-photon.pdf");
    delete c;
  }
  else if(ispion==1)
  {
    c->Print("plots/Smear-pion.pdf");
    delete c;
  }

  return;
}

void draw_Smear()
{
  //TFile *f = new TFile("data/MissingRatio-histo.root");
  TFile *f = new TFile("data/AnaPHPythia-histo.root");
  TObjArray *Glist = new TObjArray();

  for(int ispion=0; ispion<2; ispion++)
  {
    cout << "\nispion " << ispion << endl;
    GenerateSmear(f, Glist, ispion);
  }

  TFile *fout = new TFile("data/Smear.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
