void draw_Eta_Phi()
{
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  const Int_t phibin[9] = {2, 2+18, 2+18*2, 2+18*3, 4+18*4, 4+18*5, 4+18*6, 4+18*6+24, 4+18*6+24*2};

  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/total.root");
  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");

  mc(0, 2,5);

  for(Int_t sec=0; sec<8; sec++)
  {
    Int_t part;
    if(sec < 4) part = 0;
    else if(sec < 6) part = 1;
    else part = 2;

    TH2 *h2_eta_phi_data = (TH2*)f_data->Get( Form("h2_photon_eta_phi_part%d",part) );
    TH2 *h2_eta_phi_sim = (TH2*)f_sim->Get( Form("h2_pion_eta_phi_part%d",part) );

    mcd(0, sec+1);
    TH1 *h_eta_data = h2_eta_phi_data->ProjectionX("_px",phibin[sec],phibin[sec+1])->Clone("h_eta_data");
    TH1 *h_eta_sim = h2_eta_phi_sim->ProjectionX("_px",phibin[sec],phibin[sec+1])->Clone("h_eta_sim");
    Double_t scale = h_eta_data->Integral(0,-1) / h_eta_sim->Integral(0,-1);
    h_eta_sim->Scale(scale);
    h_eta_data->SetTitle( Form("#eta dist., sector %d",sec) );
    aset(h_eta_data);
    style(h_eta_data);
    h_eta_data->Draw("E");
    h_eta_sim->Draw("HIST SAME");
  }

  for(Int_t part=0; part<3; part++)
  {
    TH2 *h2_eta_phi_data = (TH2*)f_data->Get( Form("h2_photon_eta_phi_part%d",part) );
    TH2 *h2_eta_phi_sim = (TH2*)f_sim->Get( Form("h2_pion_eta_phi_part%d",part) );

    mcd(0, (part+1)/2+9);
    TH1 *h_phi_data = h2_eta_phi_data->ProjectionY()->Clone("h_phi_data");
    TH1 *h_phi_sim = h2_eta_phi_sim->ProjectionY()->Clone("h_phi_sim");
    scale = h_phi_data->Integral(0,-1) / h_phi_sim->Integral(0,-1);
    h_phi_sim->Scale(scale);
    h_phi_data->SetTitle( Form("#phi dist., part %d",part) );
    if(part == 0)
      aset(h_phi_data, "","", -0.6,1.);
    else
      aset(h_phi_data, "","", 2.14,3.8);
    style(h_phi_data);
    Double_t max_phi = h_phi_data->GetMaximum();
    h_phi_data->SetMaximum(1.1*max_phi);
    if(part < 2)
      h_phi_data->Draw("E");
    else
      h_phi_data->Draw("E SAME");
    h_phi_sim->Draw("HIST SAME");
  }

  c0->Print("plots/DirphEtaPhi.pdf");
}
