void draw_pTWeight()
{
  const Int_t npT = 30;
  Double_t apT[npT+1] = { 0., 
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };
  Double_t pT[npT];
  for(Int_t i=0; i<npT; i++)
    pT[i] = (apT[i] + apT[i+1]) / 2.;
  Double_t wpT_trk[npT] = {0,
    0, 0.307319, 0.186529, 0.137168, 0.102891, 0.0784738, 0.0574793, 0.0402374, 0.0271419, 0.0184652,
    0.0127012, 0.00886497, 0.00609304, 0.00424588, 0.00307189, 0.00218866, 0.00163727, 0.00120101, 0.000907576, 0.000483097,
    0.000187945, 8.33357e-05, 4.30282e-05, 2.11378e-05, 1.13421e-05, 6.3131e-06, 4.1332e-06, 1.9904e-06, 1.10784e-06};
  Double_t wpT_ev[npT] = {0,
    0, 0.412878, 0.208844, 0.127911, 0.0826409, 0.0554566, 0.0370308, 0.0245368, 0.0159906, 0.0105124,
    0.00711123, 0.00489351, 0.00333172, 0.00229371, 0.00162967, 0.00117611, 0.000861856, 0.000637708, 0.000480928, 0.000257418,
    9.95953e-05, 4.28212e-05, 2.10317e-05, 1.05211e-05, 5.57898e-06, 3.38423e-06, 1.81053e-06, 1.23685e-06, 7.47373e-07};
  //const Double_t gy[30] = {3.94174e+09,1e9,1e9,1e9,1e9,109923,25315.3,39417,25698,12975,7678.09,4793.76,2506.81,1537.86,960.177,596.395,399.283,251.801,192.723,132.23,61.5391,21.5451,8.59671,3.96389,2.03649,1.05427,0.560866,0.369181,0.141591,1e9,};
  //const Double_t egy[30] = {9.93204e+10,1e9,1e9,1e9,1e9,574931,4154.64,6054.26,3950,2004.45,1184.05,737.111,385.065,237.484,149.483,91.5576,61.3683,38.636,29.6415,20.3624,9.28325,3.25007,1.29766,0.598032,0.30808,0.160228,1e9,1e9,1e9,1e9,};
  const Double_t gy[30] = {4.53486e+10,5.88639e+08,1e9,793699,834434,271543,100834,54796.3,29149.8,13930.7,7795.81,4858.02,2538.79,1532.88,968.076,610.588,408.872,256.572,194.607,135.07,62.6592,21.769,9.02713,4.11342,2.08879,1.10271,0.582433,0.37722,0.162604,1e9};
  const Double_t egy[30] = {2.75178e+12,2.1439e+09,1e9,140225,133855,42870.4,15559.4,8463.26,4495.57,2140.04,1201.47,749.4,391.079,238.199,151.635,93.7945,63.0925,39.54,30.1518,20.9568,9.52342,3.31564,1.37202,0.626348,0.318893,0.170641,1e9,1e9,1e9,1e9};

  Double_t trigE[npT] = {0,0.00231696,0.00524859,0.0123286,0.0305177,0.0810258,0.17948,0.307198,0.436089,0.548318,0.637737,0.699696,0.742986,0.779692,0.801889,0.811844,0.826064,0.833043,0.837983,0.835118,0.848851,0.841841,0.821763,0.811057,0.788889,0.744472,0.675676,0.698795,0.807018,0.621622};
  Double_t acc[npT] = {0,0,0.0609309,0.0865784,0.102778,0.113609,0.11807,0.120517,0.129105,0.13699,0.132367,0.13558,0.131447,0.14864,0.143204,0.144058,0.15351,0.15129,0.15263,0.155227,0.15733,0.168737,0.171907,0.175655,0.179981,0.185073,0.185789,0.187877,0.191139,0.172771};

  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ertb/total.root");
  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");
  TH1 *h_pt = hn_2photon->Projection(2);

  //TFile *f = new TFile("AnaPHPythia-histo.root");
  //TH2 *h2_pion = (TH2*)f->Get("h2_measured_pion");
  //h2_pion->GetYaxis()->SetRange(24,24);
  //TH1 *h_pt = h2_pion->ProjectionX();

  //TFile *f = new TFile("AnaPHPythia-histo-noweight.root");
  //THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  //TH1 *h_pt = hn_pion->Projection(0);

  Double_t wpT_data[npT] = {};
  for(Int_t i=0; i<npT; i++)
  {
    wpT_data[i] = h_pt->GetBinContent(i+1);
    if(i >= 20)
      wpT_data[i] /= 4.;
    //Double_t corr = (trigE[i] * acc[i]);
    //if(corr > 0.)
    //  wpT_data[i] /= corr;
  }

  TF1 *cross1 = new TF1("cross1","372.8*x*pow(1.2282/(x+1.2282),10.00)",1,30);
  TF1 *cross2 = new TF1("cross2","x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))",1,30);
  cross2->SetParameters( 2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00 );
  //cross2->SetParameters( 229.6, 1.466, 10.654, 14.43, 8.1028, 4.5, 0.114 );

  Double_t scale_trk = wpT_trk[11] / cross2->Eval(5.75);
  Double_t scale_ev = wpT_ev[11] / cross2->Eval(5.75);
  Double_t scale_data = wpT_data[11] / cross2->Eval(5.75);
  Double_t scale_gy = gy[11] / cross2->Eval(5.75);
  for(Int_t i=0; i<npT; i++)
  {
    wpT_trk[i] /= scale_trk;
    wpT_ev[i] /= scale_ev;
    wpT_data[i] /= scale_data;
    gy[i] /= scale_gy;
    egy[i] /= scale_gy;
  }
  TGraph *gr_trk = new TGraph(npT, pT, wpT_trk);
  TGraph *gr_ev = new TGraph(npT, pT, wpT_ev);
  TGraph *gr_data = new TGraph(npT, pT, wpT_data);
  TGraphErrors *gr_gy = new TGraphErrors(npT, pT, gy, 0, egy);

  TCanvas *c = new TCanvas("c", "", 600, 600);
  gPad->SetLogy();

  cross1->SetLineColor(kRed);
  //cross1->Draw();
  cross2->SetLineColor(kBlack);
  cross2->Draw();
  gr_trk->SetMarkerColor(kBlue);
  gr_trk->SetMarkerStyle(20);
  //gr_trk->Draw("PSAME");
  gr_ev->SetMarkerColor(kYellow);
  gr_ev->SetMarkerStyle(21);
  //gr_ev->Draw("PSAME");
  gr_data->SetMarkerColor(kGreen);
  gr_data->SetMarkerStyle(22);
  //gr_data->Draw("PSAME");
  gr_gy->SetMarkerColor(kBlue);
  gr_gy->SetMarkerStyle(23);
  gr_gy->Draw("PSAME");

  c->Print("pTWeight-ertc.pdf");
}
