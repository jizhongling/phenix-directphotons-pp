TGraph** CreateGraph(TFile *f, Int_t arm)
{
  const Double_t PI = 3.1415927;
  const Double_t BBCCount = 4.193255 * pow(10,12);  // 4193255288313
  const Double_t BBCRatio = 0.2;
  const Double_t BBCCross = 32.5 * pow(10,9);
  const Double_t Acceptance = 0.201318;
  const Double_t BR[30] = {0.0444823,0.155126,0.189877,0.209775,0.227416,0.235346,0.242619,0.258304,0.254693,0.219828,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313};

  const Int_t gn = 30;
  Double_t *gx = new Double_t[gn];
  Double_t *gy[3];
  for(Int_t i=0; i<3; i++)
  {
    gy[i] = new Double_t[gn];
    for(Int_t j=0; j<gn; j++)
    {
      gx[j] = 0.;
      gy[i][j] = 0.;
    }
  }

  if(arm == 0)
  {
    const Int_t sector_low = 1;
    const Int_t sector_high = 4;
    //const Double_t MissR[30] = {1.17052, 1.02665, 1.06256, 0.869176, 0.860046, 0.771741, 0.742194, 0.687999, 0.632733, 0.792072, 0.813304, 0.699638, 0.774388, 0.657568, 0.558923, 0.732799, 0.649979, 0.54994, 0.572371, 0.476033, 0.524398, 0.658759, 0.725874, 0.477244, 0.837308, 1.13577, 1.1215, 1.69883, 1.80851, 1.98848};
    const Double_t MissR[30] = {0.576394, 1.5215, 2.09614, 1.55921, 1.98752, 1.47444, 1.54181, 1.73786, 1.36634, 1.40074, 1.27356, 1.14209, 1.11697, 1.03645, 0.886776, 1.00954, 0.928018, 0.729948, 0.753175, 0.634566, 0.928788, 1.14879, 1.18122, 0.905212, 1.30759, 1.91196, 2.2468, 3.98329, 5.55577, 7.05814};
    const Double_t TrigE[30] = {0.000164287,0.000226765,0.000554624,0.00132205,0.00273391,0.00653341,0.0205823,0.0836921,0.247432,0.506061,0.711566,0.897269,0.920601,0.926923,0.952096,0.963303,0.897059,0.916667,0.931034,1,0.930233,0.916667,1,1,1,1,1,1,1,1};
  }
  else
  {
    const Int_t sector_low = 5;
    const Int_t sector_high = 8;
    //const Double_t MissR[30] = {0.968893, 0.814013, 0.887797, 0.794183, 0.61442, 0.658933, 0.691201, 0.627664, 0.553573, 0.585781, 0.526425, 0.565531, 0.375501, 0.469886, 0.398121, 0.398228, 0.281594, 0.346736, 0.348726, 0.358025, 0.375395, 0.394456, 0.325099, 0.34044, 0.398916, 0.383546, 0.586689, 0.94373, 0.874969, 1.76884};
    const Double_t MissR[30] = {0.761121, 1.21671, 1.77948, 1.53485, 1.23296, 1.18563, 1.27982, 1.4427, 1.188, 1.22016, 0.984438, 0.974272, 0.872964, 0.72481, 0.726905, 0.57041, 0.421232, 0.448614, 0.453693, 0.491227, 0.692569, 0.688636, 0.597951, 0.62367, 0.800111, 0.857026, 1.14432, 1.86413, 2.59804, 6.26555};
    const Double_t TrigE[30] = {0.000175027,0.000228947,0.000544174,0.00116377,0.002314,0.00449421,0.0283686,0.0850631,0.195189,0.340642,0.491794,0.645161,0.753769,0.774536,0.794643,0.833333,0.837607,0.866667,0.897959,0.840909,0.860759,0.916667,1,1,1,1,1,1,1,1};
  }

  TH1 *h_events = (TH1*)f->Get("h_events");
  Double_t nevents = h_events->GetBinContent(1);
  cout << "nevents=" << nevents << endl;

  TH2 *h2_1photon = (TH2*)f->Get("h2_1photon");
  TH1 *h_1photon = h2_1photon->ProjectionY("h_1photon", sector_low, sector_high);

  TH2 *h2_sig_extra = (TH2*)f->Get("h2_sig_extra");
  TH1 *h_sig_extra = h2_sig_extra->ProjectionY("h_sig_extra", sector_low, sector_high);

  TH2 *h2_bg_extra = (TH2*)f->Get("h2_bg_extra");
  TH1 *h_bg_extra = h2_bg_extra->ProjectionY("h_bg_extra", sector_low, sector_high);

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");

  TAxis *axis0 = hn_2photon->GetAxis(0);
  TAxis *axis1 = hn_2photon->GetAxis(1);
  TAxis *axis2 = hn_2photon->GetAxis(2);

  axis0->SetRange(sector_low, sector_high);

  Int_t bin047 = axis2->FindBin(0.047);
  Int_t bin067 = axis2->FindBin(0.067);
  Int_t bin087 = axis2->FindBin(0.087);
  Int_t bin097 = axis2->FindBin(0.097);
  Int_t bin112 = axis2->FindBin(0.112);
  Int_t bin162 = axis2->FindBin(0.162);
  Int_t bin177 = axis2->FindBin(0.177);
  Int_t bin187 = axis2->FindBin(0.187);
  Int_t bin212 = axis2->FindBin(0.212);
  Int_t bin227 = axis2->FindBin(0.227);

  const Int_t nData = 256;
  vector<Double_t> x(nData), y(nData), sigma_y(nData);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c->Divide(6,5);

  Int_t ipad = 1;
  for(Int_t ipt=1; ipt<30; ipt++)
  {
    c->cd(ipad++);

    char buf[100];
    Double_t low = axis1->GetBinLowEdge(ipt+1);
    Double_t high = axis1->GetBinLowEdge(ipt+2);
    sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

    axis1->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(2);
    h_inv_mass->SetTitle(buf);

    Double_t nphoton = h_1photon->GetBinContent(ipt+1);
    Double_t nsig_extra = h_sig_extra->GetBinContent(ipt+1);
    Double_t nbg_extra = h_bg_extra->GetBinContent(ipt+1);
    Double_t nextra = nsig_extra - nbg_extra/2.;
    Double_t npair = 0.;
    Double_t nbgfit = 0.;
    Double_t nbgside = 0.;
    Double_t nbggpr = 0.;
    Double_t dnbggpr = 0.;

    x.clear();
    y.clear();
    sigma_y.clear();

    for(Int_t ib=bin067; ib<bin087; ib++)
    {
      Double_t xx = axis2->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    for(Int_t ib=bin187; ib<bin212; ib++)
    {
      Double_t xx = axis2->GetBinCenter(ib);
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

    Double_t par[10];
    h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
    fn2->SetParameters( fn1->GetParameters() );
    h_inv_mass->Fit(fn2, "Q", "", 0.047, 0.227);
    fn2->GetParameters(par);
    fn3->SetParameters(par[3], par[4], par[5], par[6]);

    for(Int_t ib=bin047; ib<bin097; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    for(Int_t ib=bin177; ib<bin227; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    nbgside /= 2.;

    for(Int_t ib=bin112; ib<bin162; ib++)
    {
      npair += h_inv_mass->GetBinContent(ib);
      Double_t bincenter = axis2->GetBinCenter(ib);
      nbgfit += fn3->Eval(bincenter);
    }

    Double_t nfit = npair - nbgfit;
    Double_t nsub = npair - nbgside;
    Double_t ngpr = npair - nbggpr;

    Double_t npion = nfit * ( 1. + MissR[ipt] );
    Double_t ndecay = npion * ( 1. + BR[ipt] );
    Double_t ndirect = nphoton - ndecay;
    gx[ipt] = (low + high) / 2.;
    gy[0][ipt] = ndirect / ndecay;
    gy[1][ipt] = BBCCross * (ndirect/BBCCount) / (2*PI*gx[ipt]) / ((high-low)*0.7) / Acceptance / TrigE[ipt] / BBCRatio;
    gy[2][ipt] = BBCCross * (npion/2/BBCCount) / (2*PI*gx[ipt]) / ((high-low)*0.7) / Acceptance / TrigE[ipt] / BBCRatio;

    cout << "pT=" << low << "-" << high << "\tnphoton=" << nphoton
      << "\tnpair=" << npair << "\tnfit=" << nfit << "\tnsub=" << nsub
      << "\tngpr=" << ngpr << endl;
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

  TGraph **graph = new TGraph*[3];
  for(Int_t i=0; i<3; i++)
    graph[i] = new TGraph(gn, gx, gy[i]);
  return graph;
}

void draw_DirectPhoton()
{
  gSystem->Load("libGausProc.so");
  gROOT->ProcessLine(".L BgGPR.C");

  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/total.root");

  TCanvas *c0 = new TCanvas("c0", "Canvas", 1800, 600);
  gStyle->SetOptStat(0);
  c0->Divide(3,1);

  TGraph **gr[2];
  for(Int_t arm=0; arm<2; arm++)
  {
    cout << "arm " << arm << endl;
    gr[arm] = CreateGraph(f, arm);
    gr[arm][0]->SetTitle("DirectPhoton/DecayPhoton");
    gr[arm][1]->SetTitle("DirectPhoton Cross Section");
    gr[arm][2]->SetTitle("#pi^{0} Cross Section");
    for(Int_t i=0; i<3; i++)
    {
      c0->cd(i+1);
      if(i==0)
      {
        gr[arm][i]->GetYaxis()->SetTitle("Ratio");
        gr[arm][i]->GetYaxis()->SetRangeUser(0., 3.);
      }
      else
      {
        gPad->SetLogy();
        gr[arm][i]->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3} [pb*GeV^{-2}*c^{-3}]");
        gr[arm][i]->GetYaxis()->SetRangeUser(0.1, pow(10,5));
      }
      gr[arm][i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      gr[arm][i]->GetYaxis()->SetTitleOffset(1.2);
      gr[arm][i]->GetXaxis()->SetRangeUser(0., 30.);
      gr[arm][i]->SetMarkerColor(arm+1);
      gr[arm][i]->SetMarkerStyle(arm+20);
      gr[arm][i]->SetMarkerSize(1.);
      if(arm == 0)
        gr[arm][i]->Draw("AP");
      else
        gr[arm][i]->Draw("P");
    }
  }

  c0->cd(1);
  TLegend *leg = new TLegend(0.1, 0.7, 0.4, 0.9);
  leg->AddEntry(gr[0][0], "West arm", "P");
  leg->AddEntry(gr[1][0], "East arm", "P");
  leg->Draw();

  c0->Print("DirectPhoton.pdf");
}
