void draw_Eta_Phi()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  const int phibin[9] = {2, 2+19, 2+19*2, 2+19*3, 4+19*4, 4+19*5, 4+19*6, 4+19*6+25, 4+19*6+25*2};

  TFile *f_data = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");
  TFile *f_sim = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo.root");

  mc(0, 4,5);

  TH2 *h2_eta_phi_data[3];
  TH2 *h2_eta_phi_sim[3];
  for(int part=0; part<3; part++)
  {
    h2_eta_phi_sim[part] = (TH2*)f_sim->Get( Form("h2_photon_eta_phi_part%d",part) );
    h2_eta_phi_data[part] = (TH2*)f_data->Get( Form("h2_eta_phi_%d",part+3) );
    //TH2 *h2_tmp = (TH2*)f_data->Get( Form("h2_eta_phi_%d",part+3) );
    //h2_eta_phi_data[part]->Add(h2_tmp);
    //delete h2_tmp;

    for(int binx=1; binx<=h2_eta_phi_data[part]->GetNbinsX(); binx++)
      for(int biny=1; biny<=h2_eta_phi_data[part]->GetNbinsY(); biny++)
      {
        if( h2_eta_phi_data[part]->GetBinContent(binx,biny) <= 0. )
          h2_eta_phi_sim[part]->SetBinContent(binx,biny,0.);
        //if( biny >= 130 && biny <= 132 )
        //{
        //  h2_eta_phi_data[part]->SetBinContent(binx,biny,0.);
        //  h2_eta_phi_sim[part]->SetBinContent(binx,biny,0.);
        //}
      }
  }

  for(int sec=0; sec<8; sec++)
  {
    int part;
    if(sec < 4) part = 0;
    else if(sec < 6) part = 1;
    else part = 2;

    mcd(0, sec+1);
    TH1 *h_eta_data =(TH1*)h2_eta_phi_data[part]->ProjectionX("_px",phibin[sec],phibin[sec+1])->Clone("h_eta_data");
    TH1 *h_eta_sim = (TH1*)h2_eta_phi_sim[part]->ProjectionX("_px",phibin[sec],phibin[sec+1])->Clone("h_eta_sim");
    double scale = h_eta_data->Integral(0,-1) / h_eta_sim->Integral(0,-1);
    if(part==2) scale *= 0.011/0.008;
    h_eta_sim->Scale(scale);
    h_eta_data->SetTitle( Form("#eta dist., sector %d",sec) );
    aset(h_eta_data);
    style(h_eta_data);
    double max_eta = h_eta_sim->GetMaximum();
    h_eta_data->SetMaximum(1.1*max_eta);
    h_eta_data->Draw("E");
    h_eta_sim->Draw("HIST SAME");
  }

  for(int part=2; part>=0; part--)
  {
    mcd(0, part+9);
    TH1 *h_eta_data =(TH1*)h2_eta_phi_data[part]->ProjectionX()->Clone("h_eta_data");
    TH1 *h_eta_sim = (TH1*)h2_eta_phi_sim[part]->ProjectionX()->Clone("h_eta_sim");
    scale = h_eta_data->Integral(0,-1) / h_eta_sim->Integral(0,-1);
    h_eta_sim->Scale(scale);
    h_eta_data->SetTitle( Form("#eta dist., part %d",part) );
    style(h_eta_data);
    double max_eta = h_eta_data->GetMaximum();
    h_eta_data->SetMaximum(1.1*max_eta);
    h_eta_data->Draw("E");
    h_eta_sim->Draw("HIST SAME");

    mcd(0, (part+1)/2+13);
    TH1 *h_phi_data =(TH1*)h2_eta_phi_data[part]->ProjectionY()->Clone("h_phi_data");
    TH1 *h_phi_sim = (TH1*)h2_eta_phi_sim[part]->ProjectionY()->Clone("h_phi_sim");
    double scale = h_phi_data->Integral(0,-1) / h_phi_sim->Integral(0,-1);
    h_phi_sim->Scale(scale);
    if(part == 2)
    {
      double sc2gl = 0.011 / 0.008;
      h_phi_data->Scale(sc2gl);
      h_phi_sim->Scale(sc2gl);
    }
    //for(int binx=1; binx<=h_phi_data->GetNbinsX(); binx++)
    //  if( h_phi_data->GetBinContent(binx) > 15000. )
    //  {
    //    h_phi_data->SetBinContent(binx,0.);
    //    h_phi_sim->SetBinContent(binx,0.);
    //  }
    h_phi_data->SetTitle( Form("#phi dist., part %d",part) );
    if(part == 0)
      aset(h_phi_data, "","", -0.6,1.);
    else
      aset(h_phi_data, "","", 2.14,3.8);
    style(h_phi_data);
    double max_phi = h_phi_data->GetMaximum();
    h_phi_data->SetMaximum(1.1*max_phi);
    if(part != 1)
      h_phi_data->Draw("E");
    else
      h_phi_data->Draw("E SAME");
    h_phi_sim->Draw("HIST SAME");

    mcd(0, part*2+15);
    h2_eta_phi_data[part]->SetTitle( Form("#eta and #phi from data for part %d",part) );
    h2_eta_phi_data[part]->GetYaxis()->SetRange(phibin[secl[part]-1], phibin[sech[part]]-3);
    h2_eta_phi_data[part]->Draw("COLZ");

    mcd(0, part*2+16);
    h2_eta_phi_sim[part]->SetTitle( Form("#eta and #phi from FastMC for part %d",part) );
    h2_eta_phi_sim[part]->GetYaxis()->SetRange(phibin[secl[part]-1], phibin[sech[part]]-3);
    h2_eta_phi_sim[part]->Draw("COLZ");
  }

  c0->Print("plots/DirphEtaPhi-isolated.pdf");
}
