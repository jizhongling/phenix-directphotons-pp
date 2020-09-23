void draw_DCZedPhi()
{
  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-DC3sigma.root");
  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/HadronResponse-histo-minbias.root");

  TH3 *h3_dclive_data = (TH3*)f_data->Get("h3_dclive_1");
  //h3_dclive_data->Add( (TH3*)f_data->Get("h3_dclive_0") );
  h3_dclive_data->GetZaxis()->SetRange(8,30);  // pT
  TH2 *h2_phized_data = (TH2*)h3_dclive_data->Project3D("yx")->Clone("h2_phized_data");
  TH1 *h_zed_data = (TH1*)h2_phized_data->ProjectionX("h_zed_data", 0,-1)->Clone("h_zed_data");
  TH1 *h_phi_data = (TH1*)h2_phized_data->ProjectionY("h_phi_data", 0,-1)->Clone("h_phi_data");;

  THnSparse *hn_dclive_sim = (THnSparse*)f_sim->Get("hn_dclive");
  hn_dclive_sim->GetAxis(4)->SetRange(2,2);  // ERT
  hn_dclive_sim->GetAxis(3)->SetRange(2,2);  // isDCGood
  hn_dclive_sim->GetAxis(2)->SetRange(8,30);  // pT
  TH2 *h2_phized_sim = hn_dclive_sim->Projection(1,0);

  //mc();
  //mcd();
  //const int sec[] = {1, 25, 40, 50};
  //for(int part=0; part<3; part++)
  //{
  //  hn_dclive_sim->GetAxis(1)->SetRange(sec[part],sec[part+1]);  // Sector
  //  hn_dclive_sim->GetAxis(3)->SetRange(2,2);  // isDCGood
  //  hn_dclive_sim->GetAxis(4)->SetRange(1,2);  // ERT
  //  TH1 *h_total = hn_dclive_sim->Projection(2);
  //  hn_dclive_sim->GetAxis(4)->SetRange(2,2);  // ERT
  //  TH1 *h_passed = hn_dclive_sim->Projection(2);
  //  TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_passed, h_total, "n");
  //  aset(gr, "p_{T} [GeV/c]","", 3.1,15., 0.,0.2);
  //  style(gr, part+20, part+1);
  //  if(part==0)
  //    gr->Draw("APE");
  //  else
  //    gr->Draw("PE");
  //}
  //return;

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
