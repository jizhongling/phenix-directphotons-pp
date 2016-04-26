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

    //Double_t nphoton = h_1photon->GetBinContent(ipt+1);
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
    TF1 *fn2 = new TF1("fn2", "gaus(0)+pol2(3)", 0., 0.5);
    TF1 *fn3 = new TF1("fn3", "pol2", 0., 0.5);

    Double_t par[6];
    h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
    fn2->SetParameters( fn1->GetParameters() );
    h_inv_mass->Fit(fn2, "Q0", "", 0.047, 0.227);
    fn2->GetParameters(par);
    fn3->SetParameters(par[3], par[4], par[5]);

    for(Int_t ib=bin047; ib<bin097; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    for(Int_t ib=bin177; ib<bin227; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    nbgside /= 2.;

    for(Int_t ib=bin112; ib<bin162; ib++)
    {
      npair += h_inv_mass->GetBinContent(ib);
      Double_t bincenter = axis1->GetBinCenter(ib);
      nbgfit += fn3->Eval(bincenter);
    }

    Double_t nfit = npair - nbgfit;
    Double_t nsub = npair - nbgside;
    Double_t ngpr = npair - nbggpr;

    cout << "pT=" << low << "-" << high << "\tnpair=" << npair
      << "\tnfit=" << nfit << "\tnsub=" << nsub
      << "\tngpr=" << ngpr << "\tdnbggpr=" << dnbggpr << endl;

    h_inv_mass->Delete();
  }

  //c->Print("DirectPhoton.pdf");
}
