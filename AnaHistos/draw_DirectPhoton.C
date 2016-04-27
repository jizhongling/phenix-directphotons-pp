void draw_arm(TFile *f, int arm)
{
  if(arm == 0)
  {
    const int sector_low = 1;
    const int sector_high = 4;
    const double MissR[30] = {1.17052, 1.02665, 1.06256, 0.869176, 0.860046, 0.771741, 0.742194, 0.687999, 0.632733, 0.792072, 0.813304, 0.699638, 0.774388, 0.657568, 0.558923, 0.732799, 0.649979, 0.54994, 0.572371, 0.476033, 0.524398, 0.658759, 0.725874, 0.477244, 0.837308, 1.13577, 1.1215, 1.69883, 1.80851, 1.98848};
  }
  else
  {
    const int sector_low = 5;
    const int sector_high = 8;
    const double MissR[30] = {0.968893, 0.814013, 0.887797, 0.794183, 0.61442, 0.658933, 0.691201, 0.627664, 0.553573, 0.585781, 0.526425, 0.565531, 0.375501, 0.469886, 0.398121, 0.398228, 0.281594, 0.346736, 0.348726, 0.358025, 0.375395, 0.394456, 0.325099, 0.34044, 0.398916, 0.383546, 0.586689, 0.94373, 0.874969, 1.76884};
  }

  TH1 *h_events = (TH1*)f->Get("event_counter");
  Double_t nevents = h_events->GetBinContent(4);
  cout << "nevents=" << nevents << endl;

  THnSparse *hn_1photon = (THnSparse*)f->Get("pt_1photon");
  hn_1photon->GetAxis(2)->SetRange(1,1);
  hn_1photon->GetAxis(1)->SetRange(sector_low, sector_high);
  TH1 *h_1photon = hn_1photon->Projection(0);

  THnSparse *hn_2photon = (THnSparse*)f->Get("inv_mass_2photon");

  TAxis *axis0 = hn_2photon->GetAxis(0);
  TAxis *axis1 = hn_2photon->GetAxis(1);
  TAxis *axis2 = hn_2photon->GetAxis(2);
  TAxis *axis3 = hn_2photon->GetAxis(3);

  axis3->SetRange(2,2);
  axis2->SetRange(sector_low, sector_high);

  Int_t bin047 = axis1->FindBin(0.047);
  Int_t bin067 = axis1->FindBin(0.067);
  Int_t bin087 = axis1->FindBin(0.087);
  Int_t bin097 = axis1->FindBin(0.097);
  Int_t bin112 = axis1->FindBin(0.112);
  Int_t bin162 = axis1->FindBin(0.162);
  Int_t bin177 = axis1->FindBin(0.177);
  Int_t bin187 = axis1->FindBin(0.187);
  Int_t bin212 = axis1->FindBin(0.212);
  Int_t bin227 = axis1->FindBin(0.227);

  const Int_t nData = 256;
  vector<Double_t> x(nData), y(nData), sigma_y(nData);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2400);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c->Divide(4,4);

  Int_t ipad = 1;
  for(Int_t ipt=5; ipt<21; ipt++)
  {
    c->cd(ipad++);

    char buf[100];
    Double_t low = axis0->GetBinLowEdge(ipt+1);
    Double_t high = axis0->GetBinLowEdge(ipt+2);
    sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

    axis0->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(1);
    h_inv_mass->SetTitle(buf);
    h_inv_mass->DrawCopy();

    Double_t nphoton = h_1photon->GetBinContent(ipt+1);
    Double_t npair = 0.;
    Double_t nbgfitgaus = 0.;
    Double_t nbgfitvoigt = 0.;
    Double_t nbgside = 0.;
    Double_t nbggpr = 0.;
    Double_t dnbggpr = 0.;

    x.clear();
    y.clear();
    sigma_y.clear();

    for(Int_t ib=bin067; ib<bin087; ib++)
    {
      Double_t xx = axis1->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    for(Int_t ib=bin187; ib<bin212; ib++)
    {
      Double_t xx = axis1->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    BgGPR(x, y, sigma_y, nbggpr, dnbggpr);
    nbggpr /= 0.001;
    dnbggpr /= 0.001;

    TF1 *fn1 = new TF1("fn1", "gaus", 0., 0.5);
    TF1 *fn2 = new TF1("fn2", "gaus(0)+pol3(3)", 0., 0.5);
    TF1 *fn3 = new TF1("fn3", "pol3", 0., 0.5);

    TF1 *fn4 = new TF1("fn4", "[0]*TMath::Voigt(x-[1],[2],[3])", 0., 0.5);
    TF1 *fn5 = new TF1("fn5", "[0]*TMath::Voigt(x-[1],[2],[3])+pol3(4)", 0., 0.5);
    TF1 *fn6 = new TF1("fn6", "pol3", 0., 0.5);

    Double_t par1[10];
    Double_t par2[10];

    h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
    fn2->SetParameters( fn1->GetParameters() );
    h_inv_mass->Fit(fn2, "Q0", "", 0.047, 0.227);
    fn2->SetLineStyle(1);
    fn2->SetLineColor(kRed);
    fn2->SetRange(0.047,0.227);
    fn2->Draw("CSAME");
    fn2->GetParameters(par1);
    fn3->SetParameters(par1[3], par1[4], par1[5], par1[6]);

    fn4->SetParameters(1000., 0.137, 0.005, 0.005);
    h_inv_mass->Fit(fn4, "Q0", "", 0.112, 0.162);
    fn5->SetParameters( fn4->GetParameters() );
    h_inv_mass->Fit(fn5, "Q0", "", 0.047, 0.227);
    fn5->SetLineStyle(3);
    fn5->SetLineColor(kGreen);
    fn5->SetRange(0.047,0.227);
    fn5->Draw("CSAME");
    fn5->GetParameters(par2);
    fn6->SetParameters(par2[4], par2[5], par2[6], par2[7]);

    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(fn2, "gaus+pol3", "L");
    leg->AddEntry(fn5, "voigt+pol3", "L");
    leg->Draw();

    for(Int_t ib=bin047; ib<bin097; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    for(Int_t ib=bin177; ib<bin227; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    nbgside /= 2.;

    for(Int_t ib=bin112; ib<bin162; ib++)
    {
      npair += h_inv_mass->GetBinContent(ib);
      Double_t bincenter = axis1->GetBinCenter(ib);
      nbgfitgaus += fn3->Eval(bincenter);
      nbgfitvoigt += fn6->Eval(bincenter);
    }

    Double_t nfitgaus = npair - nbgfitgaus;
    Double_t nfitvoigt = npair - nbgfitvoigt;
    Double_t nsub = npair - nbgside;
    Double_t ngpr = npair - nbggpr;

    cout << "pT=" << low << "-" << high << "\tnphoton=" << nphoton
      << "\tnpair=" << npair << "\tnfitgaus=" << nfitgaus << "\tnfitvoigt=" << nfitvoigt
      << "\tnsub=" << nsub << "\tngpr=" << ngpr << endl;

    h_inv_mass->Delete();
  }

  if(arm == 0)
  {
    c->Print("DirectPhoton-west.pdf");
    delete c;
  }
  else
  {
    c->Print("DirectPhoton-east.pdf");
    delete c;
  }
}

void draw_DirectPhoton()
{
  gSystem->Load("libGausProc.so");
  gROOT->ProcessLine(".L BgGPR.C");

  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/8511/data/total.root");

  cout << "West arm" << endl;
  draw_arm(f, 0);

  cout << "East arm" << endl;
  draw_arm(f, 1);
}
