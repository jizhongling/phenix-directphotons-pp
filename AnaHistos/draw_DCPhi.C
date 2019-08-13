void draw_DCPhi()
{
  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");
  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/HadronResponse-histo.root");

  TH3 *h3_dclive_data = (TH3*)f_data->Get("h3_dclive_1");
  TH1 *h_phi_data = h3_dclive_data->ProjectionY("h_phi_data", 0,-1, 7,30);

  THnSparse *hn_dclive_sim = (THnSparse*)f_sim->Get("hn_dclive");
  hn_dclive_sim->GetAxis(3)->SetRange(2,2);
  hn_dclive_sim->GetAxis(2)->SetRange(7,30);
  TH1 *h_phi_sim = hn_dclive_sim->Projection(1);

  double scale = h_phi_data->Integral(1,50) / h_phi_sim->Integral(1,50);
  h_phi_sim->Scale(scale);

  mc(0);
  mcd(0);
  aset(h_phi_data);
  style(h_phi_data, 20, 2);
  style(h_phi_sim, 20, 1);
  h_phi_data->Draw("HISTO");
  h_phi_sim->Draw("HISTO SAME");
  c0->Print("plots/DCPhi.pdf");
}
