int draw_photons()
{

  /* Plots to make at sqrt(s) = 200 GeV and sqrt(s) = 510 GeV:
   * pT_gamma spectrum for Direct Photons and Pi0 Decay Photons and All Photons
   * Isolated / all ratio vs pT_gamma for Direct Photons and Pi0 Decay Photons
   * Isolated / all ratio vs pT_gamma for pi0 using different cuts (energy only, energy + momentum )
   * Isolated / all ratio vs pT_gamma All photons, Direct photon, fragment photon (parent is quark)
   */

  TFile *fin_all = TFile::Open("anaphpythia_all.root");

  /* Define cuts */
  {
    TCut cut_acceptance("abs(eta) < 0.35 && pt > 1");

    TCut cut_isolation_c1("etot < 3");
  }

  /* Define histograms */
  {
    TH1F* h_pt_all = new TH1F("h_pt_all", ";p_{T}^{#gamma} (GeV/c); # entries", 25, 0.5, 25.5);
    TH1F* h_pt_pi0 = (TH1F*)h_pt_all->Clone("h_pt_pi0");
    TH1F* h_pt_dir = (TH1F*)h_pt_all->Clone("h_pt_dir");

    h_pt_all->Sumw2();
    h_pt_pi0->Sumw2();
    h_pt_dir->Sumw2();

    h_pt_all->SetLineColor(kBlue);
    h_pt_pi0->SetLineColor(kRed);
    h_pt_dir->SetLineColor(kGreen);

    TH1F* h_pt_all_pass_c1 = (TH1F*)h_pt_all->Clone("h_pt_all_pass_c1");
    TH1F* h_pt_pi0_pass_c1 = (TH1F*)h_pt_pi0->Clone("h_pt_pi0_pass_c1");
    TH1F* h_pt_dir_pass_c1 = (TH1F*)h_pt_dir->Clone("h_pt_dir_pass_c1");
  }

  TCanvas *ctemp = new TCanvas();

  /* PLOT: pT_gamma spectrum for Direct Photons and Pi0 Decay Photons and All Photons */
  {
    ctemp->cd();

    TLegend* legend_c1 = new TLegend(0.5,0.7,0.9,0.9);
    legend_c1->AddEntry(h_pt_all,"All","l");
    legend_c1->AddEntry(h_pt_pi0,"Pi0","l");
    legend_c1->AddEntry(h_pt_dir,"Direct photon","lep");

    // all
    event_truth->Draw("pt >> h_pt_all",cut_acceptance);

    // pi0's
    event_truth->Draw("pt >> h_pt_pi0",cut_acceptance && "ispromptphoton==0&&parentid==111");

    // prompt
    event_truth->Draw("pt >> h_pt_dir",cut_acceptance && "ispromptphoton==1");

    TCanvas *c1 = new TCanvas();
    c1->SetLogy();
    h_pt_all->Draw();
    h_pt_pi0->Draw("sames");
    h_pt_dir->Draw("sames");
    legend_c1->Draw();
    gPad->RedrawAxis();
  }

  /* PLOT: */
  {
    ctemp->cd();

    event_truth->Draw("pt >> h_pt_all_pass_c1",cut_acceptance && cut_isolation_c1);
    event_truth->Draw("pt >> h_pt_pi0_pass_c1",cut_acceptance && cut_isolation_c1 && "ispromptphoton==0&&parentid==111");
    event_truth->Draw("pt >> h_pt_dir_pass_c1",cut_acceptance && cut_isolation_c1 && "ispromptphoton==1");

    h_pt_all_pass_c1->Divide(h_pt_all);
    h_pt_pi0_pass_c1->Divide(h_pt_pi0);
    h_pt_dir_pass_c1->Divide(h_pt_dir);

    TCanvas *c2 = new TCanvas();
    h_pt_all_pass_c1->Draw("");
    h_pt_pi0_pass_c1->Draw("sames");
    h_pt_dir_pass_c1->Draw("sames");
    legend_c1->Draw();
    gPad->RedrawAxis();
  }

  return 0;
}
