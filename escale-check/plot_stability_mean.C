plot_stability_mean(){

  gStyle->SetOptStat(0);

  TTree *tfit = new TTree();
  tfit->ReadFile("pi0peak_fit_Run13pp510ERT_reformat.txt","run_index/F:run_number:nevents:sector:mean:dmean:sigma:dsigma:chisquare");

  TTree *tfit_raw = new TTree();
  tfit_raw->ReadFile("pi0peak_fit_raw_Run13pp510ERT_reformat.txt","run_index/F:run_number:nevents:sector:mean:dmean:sigma:dsigma:chisquare");

  TLine *lpi0 = new TLine( 0, 0.137, tfit->GetEntries("sector==0"), 0.137 );
  lpi0->SetLineColor(kRed);

  // frame
  TH1F* hframe = new TH1F("hframe","",875,0,875);
  hframe->GetYaxis()->SetRangeUser(0.12,0.16);
  hframe->GetXaxis()->SetTitle("run");
  hframe->GetYaxis()->SetTitle("mean [GeV]");
  hframe->SetLineColor(kWhite);

  // sectors PbSc-West
  tfit->Draw("mean:run_index:dmean","sector==0");
  TGraphErrors *g_mean_PbScW = new TGraphErrors( tfit->GetEntries("sector==0"), tfit->GetV2(), tfit->GetV1(), 0, tfit->GetV3());

  tfit_raw->Draw("mean:run_index:dmean","sector==0");
  TGraphErrors *g_mean_PbScW_raw = new TGraphErrors( tfit_raw->GetEntries("sector==0"), tfit_raw->GetV2(), tfit_raw->GetV1(), 0, tfit_raw->GetV3());

  // sectors PbSc-East
  tfit->Draw("mean:run_index:dmean","sector==1");
  TGraphErrors *g_mean_PbScE = new TGraphErrors( tfit->GetEntries("sector==1"), tfit->GetV2(), tfit->GetV1(), 0, tfit->GetV3());

  tfit_raw->Draw("mean:run_index:dmean","sector==1");
  TGraphErrors *g_mean_PbScE_raw = new TGraphErrors( tfit_raw->GetEntries("sector==1"), tfit_raw->GetV2(), tfit_raw->GetV1(), 0, tfit_raw->GetV3());

  // sectors PbGl-East
  tfit->Draw("mean:run_index:dmean","sector==2");
  TGraphErrors *g_mean_PbGlE = new TGraphErrors( tfit->GetEntries("sector==2"), tfit->GetV2(), tfit->GetV1(), 0, tfit->GetV3());

  tfit_raw->Draw("mean:run_index:dmean","sector==2");
  TGraphErrors *g_mean_PbGlE_raw = new TGraphErrors( tfit_raw->GetEntries("sector==2"), tfit_raw->GetV2(), tfit_raw->GetV1(), 0, tfit_raw->GetV3());


  // Plotting
  TCanvas *c_PbScW = new TCanvas();
  hframe->Draw();
  g_mean_PbScW->Draw("Psame");
  lpi0->Draw("same");
  c_PbScW->Print("plots-escale-check/mpi0_mean_PbScW.png");

  TCanvas *c_PbScW_raw = new TCanvas();
  hframe->Draw();
  g_mean_PbScW_raw->Draw("Psame");
  lpi0->Draw("same");
  c_PbScW_raw->Print("plots-escale-check/mpi0_mean_PbScW_raw.png");

  TCanvas *c_PbScE = new TCanvas();
  hframe->Draw();
  g_mean_PbScE->Draw("Psame");
  lpi0->Draw("same");
  c_PbScE->Print("plots-escale-check/mpi0_mean_PbScE.png");

  TCanvas *c_PbScE_raw = new TCanvas();
  hframe->Draw();
  g_mean_PbScE_raw->Draw("Psame");
  lpi0->Draw("same");
  c_PbScE_raw->Print("plots-escale-check/mpi0_mean_PbScE_raw.png");

  TCanvas *c_PbGlE = new TCanvas();
  hframe->Draw();
  g_mean_PbGlE->Draw("Psame");
  lpi0->Draw("same");
  c_PbGlE->Print("plots-escale-check/mpi0_mean_PbGlE.png");

  TCanvas *c_PbGlE_raw = new TCanvas();
  hframe->Draw();
  g_mean_PbGlE_raw->Draw("Psame");
  lpi0->Draw("same");
  c_PbGlE_raw->Print("plots-escale-check/mpi0_mean_PbGlE_raw.png");

}

