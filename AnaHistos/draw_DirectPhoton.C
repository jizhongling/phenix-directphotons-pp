void draw_DirectPhoton()
{
  gSystem->Load("libGausProc.so");
  gROOT->ProcessLine(".L BgGPR.C");

  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/8511/data/total.root");
  THnSparse *hn_2photon = (THnSparse*)f->Get("inv_mass_2photon");
  //TH2 *h2_1photon = (TH2*)f->Get("pT_1cluster");
  //TH1 *h_1photon = (TH1*)h2_1photon->ProjectionX("pT", 1, 4);

  TAxis *axis0 = hn_2photon->GetAxis(0);
  TAxis *axis1 = hn_2photon->GetAxis(1);
  TAxis *axis2 = hn_2photon->GetAxis(2);
  TAxis *axis3 = hn_2photon->GetAxis(3);

  axis3->SetRange(2,2);
  axis2->SetRange(1,4);

  Int_t bin047 = axis1->FindBin(0.047);
  Int_t bin097 = axis1->FindBin(0.097);
  Int_t bin112 = axis1->FindBin(0.112);
  Int_t bin162 = axis1->FindBin(0.162);
  Int_t bin177 = axis1->FindBin(0.177);
  Int_t bin227 = axis1->FindBin(0.227);

  const Int_t nData = 256;
  vector<Double_t> x1(nData), y1(nData), sigma_y1(nData);
  vector<Double_t> x2(nData), y2(nData), sigma_y2(nData);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2400);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c->Divide(4,4);

  Int_t ipad = 1;
  for(Int_t ipt=0; ipt<30; ipt++)
  {
    c->cd(ipad++);

    char buf[100];
    Double_t low = axis0->GetBinLowEdge(ipt+1);
    Double_t high = axis0->GetBinLowEdge(ipt+2);
    sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

    axis0->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = (TH1*)hn_2photon->Projection(1);
    h_inv_mass->SetTitle(buf);

    x1.clear();
    y1.clear();
    sigma_y1.clear();
    x2.clear();
    y2.clear();
    sigma_y2.clear();

    for(Int_t ib=bin047; ib<bin227; ib++)
    {
      Double_t xx = axis1->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      if( ib >= bin112 && ib < bin162 )
      {
        x1.push_back(xx);
        y1.push_back(yy);
        sigma_y1.push_back(sigma_yy);
      }
      x2.push_back(xx);
      y2.push_back(yy);
      sigma_y2.push_back(sigma_yy);
    }

    Double_t rho = BgGPR(x1, y1, sigma_y1, x2, y2, sigma_y2);

    TF1 *fn1 = new TF1("fn1", "gaus", 0., 0.5);
    TF1 *fn2 = new TF1("fn2", "gaus(0)+pol2(3)", 0., 0.5);
    TF1 *fn3 = new TF1("fn3", "pol2", 0., 0.5);

    Double_t par[6];
    h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
    fn2->SetParameters( fn1->GetParameters() );
    h_inv_mass->Fit(fn2, "QR");
    fn2->GetParameters(par);
    fn3->SetParameters(par[3], par[4], par[5]);

    //Double_t nphoton = h_1photon->GetBinContent(ipt+1);
    Double_t npair = 0.;
    Double_t nbgfit = 0.;
    Double_t nbgint = 0.;

    for(Int_t ib=bin047; ib<bin097; ib++)
      nbgint += h_inv_mass->GetBinContent(ib);
    for(Int_t ib=bin177; ib<bin227; ib++)
      nbgint += h_inv_mass->GetBinContent(ib);
    nbgint /= 2.;

    for(Int_t ib=bin112; ib<bin162; ib++)
    {
      npair += h_inv_mass->GetBinContent(ib);
      Double_t bincenter = axis1->GetBinCenter(ib);
      nbgfit += fn3->Eval(bincenter);
    }

    Double_t nfit = npair - nbgfit;
    Double_t nint = npair - nbgint;
    Double_t nintcorr = nint / (1-rho);

    cout << "pT=" << low << "-" << high << "\tnpair=" << npair
      << "\tnfit=" << nfit << "\tnint=" << nint
      << "\tnintcorr=" << nintcorr << "\trho=" << rho << endl;

    h_inv_mass->Delete();
  }

  //c->Print("DirectPhoton.pdf");
}
