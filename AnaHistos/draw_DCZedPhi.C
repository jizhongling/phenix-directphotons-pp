void draw_DCZedPhi()
{
  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");
  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/HadronResponse-histo.root");

  TH3 *h3_dclive_data = (TH3*)f_data->Get("h3_dclive_1");
  //h3_dclive_data->Add( (TH3*)f_data->Get("h3_dclive_0") );
  h3_dclive_data->GetZaxis()->SetRange(7,30);
  TH2 *h2_phized_data = (TH2*)h3_dclive_data->Project3D("yx")->Clone("h2_phized_data");
  TH1 *h_zed_data = (TH1*)h2_phized_data->ProjectionX("h_zed_data", 0,-1)->Clone("h_zed_data");
  TH1 *h_phi_data = (TH1*)h2_phized_data->ProjectionY("h_phi_data", 0,-1)->Clone("h_phi_data");;

  THnSparse *hn_dclive_sim = (THnSparse*)f_sim->Get("hn_dclive");
  hn_dclive_sim->GetAxis(3)->SetRange(2,2);
  hn_dclive_sim->GetAxis(2)->SetRange(7,30);
  TH2 *h2_phized_sim = hn_dclive_sim->Projection(1,0);
  TH1 *h_zed_sim = (TH1*)h2_phized_sim->ProjectionX("h_zed_sim", 0,-1)->Clone("h_zed_sim");
  TH1 *h_phi_sim = (TH1*)h2_phized_sim->ProjectionY("h_phi_sim", 0,-1)->Clone("h_phi_sim");;

  double scale_phized = h2_phized_data->GetEntries() / h2_phized_sim->GetEntries();
  h2_phized_sim->Scale(scale_phized);
  double scale_zed = h_zed_data->Integral(1,200) / h_zed_sim->Integral(1,200);
  h_zed_sim->Scale(scale_zed);
  double scale_phi = h_phi_data->Integral(1,25) / h_phi_sim->Integral(1,25);
  h_phi_sim->Scale(scale_phi);

  mc(0, 2,2);

  mcd(0, 1);
  aset(h_zed_sim);
  h_zed_data->SetLineColor(kRed);
  h_zed_sim->SetLineColor(kBlack);
  h_zed_sim->Draw("HISTO");
  h_zed_data->Draw("HISTO SAME");

  mcd(0, 2);
  aset(h_phi_data);
  h_phi_data->SetLineColor(kRed);
  h_phi_sim->SetLineColor(kBlack);
  h_phi_data->Draw("HISTO");
  h_phi_sim->Draw("HISTO SAME");

  mcd(0, 3);
  h2_phized_data->SetTitle("Data zed and #phi distribution");
  h2_phized_data->Draw("COLZ");

  mcd(0, 4);
  h2_phized_sim->SetTitle("PISA zed and #phi distribution");
  h2_phized_sim->Draw("COLZ");

  c0->Print("plots/DCZedPhi.pdf");
}
