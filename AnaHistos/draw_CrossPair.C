void draw_arm(TFile *f, int arm)
{
  if(arm == 0)
  {
    const int sector_low = 1;
    const int sector_high = 4;
  }
  else
  {
    const int sector_low = 5;
    const int sector_high = 8;
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
    Double_t ndiff = fabs(nbgfit - nbgside - nextra);
    Double_t ngpr = npair - nbggpr;
    Double_t nerror = fabs(nfit - ngpr);

    cout << "pT=" << low << "-" << high << "\tnphoton=" << nphoton
      << "\tnpair=" << npair << "\tnfit=" << nfit << "\tnsub=" << nsub
      << "\tndiff=" << ndiff/nfit*100. << "%\tngpr=" << ngpr
      << "\tnerror=" << nerror/nfit*100. << "%" << endl;
  }

  if(arm == 0)
  {
    c->Print("CrossPair-west.pdf");
    delete c;
  }
  else
  {
    c->Print("CrossPair-east.pdf");
    delete c;
  }
}

void draw_CrossPair()
{
  gSystem->Load("libGausProc.so");
  gROOT->ProcessLine(".L BgGPR.C");

  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/total.root");

  cout << "West arm" << endl;
  draw_arm(f, 0);

  cout << "East arm" << endl;
  draw_arm(f, 1);
}
