void test_fits()
{
  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ertb-cv/total.root");

  double sector_low = 0;
  double sector_high = 4;

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");

  /* list axis names */
  for ( unsigned axis = 0; axis < 5; axis++ )
    {
      cout << "Axis " << axis << ": " << hn_2photon->GetAxis( axis )->GetTitle() << endl;
    }

  bool ispion = true;

  TAxis *axis0 = hn_2photon->GetAxis(0);
  if(ispion == 0)
    TAxis *axis1_2 = hn_2photon->GetAxis(1);
  else if(ispion == 1)
    TAxis *axis1_2 = hn_2photon->GetAxis(2);
  TAxis *axis3 = hn_2photon->GetAxis(3);

  cout << "Set range (axis " << axis0->GetTitle() << ") to: " << sector_low << " , " << sector_high << endl;
  axis0->SetRange(sector_low, sector_high);

  int ipt = 12;
  char buf[100];
  Double_t low = axis1_2->GetBinLowEdge(ipt+1);
  Double_t high = axis1_2->GetBinLowEdge(ipt+2);
  sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

  axis1_2->SetRange(ipt+1,ipt+1);
  TH1 *h_inv_mass = hn_2photon->Projection(3);
  h_inv_mass->SetTitle(buf);

  TCanvas *c1 = new TCanvas();

  gStyle->SetOptFit(11111);

  h_inv_mass->Draw();

  TF1 *fn1 = new TF1("fn1", "gaus", 0., 0.5);
  TF1 *fn2 = new TF1("fn2", "gaus(0)+pol3(3)", 0., 0.5);
  TF1 *fn3 = new TF1("fn3", "pol3", 0., 0.5);

  Double_t par[10];
  h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
  fn2->SetParameters( fn1->GetParameters() );
  h_inv_mass->Fit(fn2, "LRSMEV", "", 0.047, 0.227);
  h_inv_mass->Fit(fn2, "LRSMEV", "", 0.047, 0.227);
  fn2->GetParameters(par);
  fn1->SetParameters(par[0], par[1], par[2]);
  fn3->SetParameters(par[3], par[4], par[5], par[6]);

  fn1->SetLineColor(kGreen);
  fn3->SetLineColor(kBlue);

  fn1->Draw("same");
  fn3->Draw("same");

}
