int draw_photons()
{

  /* Plots to make at sqrt(s) = 200 GeV and sqrt(s) = 510 GeV:
   * pT_gamma spectrum for Direct Photons and Pi0 Decay Photons and All Photons
   * Isolated / all ratio vs pT_gamma for Direct Photons and Pi0 Decay Photons
   * Isolated / all ratio vs pT_gamma for pi0 using different cuts (energy only, energy + momentum )
   * Isolated / all ratio vs pT_gamma All photons, Direct photon, fragment photon (parent is quark)
   */

  /* Open input file(s) */
  TChain* event_truth = new TChain("event_truth");
  //  event_truth->AddFile("anaphpythia_510GeV_mb.root");
  event_truth->AddFile("anaphpythia_200GeV_mb.root");

  gStyle->SetOptStat(0);

  /* Print isolation cut cone sizes */
  cout << "Isolation cut cone sizes: " << endl;
  event_truth->Scan("iso_conesize", "Entry$ ==0");

  /* Define cuts */
  {
    TCut cut_acceptance("abs(eta) < 0.35 && pt > 1");

    TCut cut_pizerophoton("ispromptphoton==0&&parentid==111");

    TCut cut_promptphoton("ispromptphoton==1");

    /* select arrray of 'isolation cut settings' */
    TObjArray arr_cut_isolation;
    arr_cut_isolation.Add( new TCut("1") );
    arr_cut_isolation.Add( new TCut("(iso_eemcal[0]+iso_ptrack[0]) < 0.1 * etot") );
    arr_cut_isolation.Add( new TCut("(iso_eemcal[1]+iso_ptrack[1]) < 0.1 * etot") );
    arr_cut_isolation.Add( new TCut("(iso_eemcal[2]+iso_ptrack[2]) < 0.1 * etot") );
    arr_cut_isolation.Add( new TCut("(iso_eemcal[3]+iso_ptrack[3]) < 0.1 * etot") );
    arr_cut_isolation.Add( new TCut("(iso_eemcal[4]+iso_ptrack[4]) < 0.1 * etot") );
  }

  /* Define and fill histograms */
  {
    TObjArray arr_hist_pt_all;
    TObjArray arr_hist_pt_pi0;
    TObjArray arr_hist_pt_dir;

    TH1F* h_pt = new TH1F("h_pt", ";p_{T}^{#gamma} (GeV/c); # entries", 25, 0.5, 25.5);
    h_pt->Sumw2();

    /* loop over isolation cut settings */
    for ( unsigned i = 0; i < arr_cut_isolation.GetEntries(); i++ )
      {
	cout << "Use Isolation Cut: " << arr_cut_isolation.At(i)->GetTitle() << endl;

        TH1F* h_pt_all = (TH1F*)h_pt->Clone( ( TString("h_pt_all_") += i ) );
        TH1F* h_pt_pi0 = (TH1F*)h_pt->Clone( ( TString("h_pt_pi0_") += i ) );
        TH1F* h_pt_dir = (TH1F*)h_pt->Clone( ( TString("h_pt_dir_") += i ) );

        h_pt_all->SetTitle( TString("Cut: ") + arr_cut_isolation.At(i)->GetTitle() );
        h_pt_pi0->SetTitle( TString("Cut: ") + arr_cut_isolation.At(i)->GetTitle() );
        h_pt_dir->SetTitle( TString("Cut: ") + arr_cut_isolation.At(i)->GetTitle() );

        h_pt_all->SetLineColor(kBlue);
        h_pt_pi0->SetLineColor(kRed);
        h_pt_dir->SetLineColor(kGreen);

        event_truth->Draw(( TString("pt >> h_pt_all_") += i ), cut_acceptance && arr_cut_isolation.At(i)->GetTitle() );
        event_truth->Draw(( TString("pt >> h_pt_pi0_") += i ), cut_acceptance && cut_pizerophoton && arr_cut_isolation.At(i)->GetTitle() );
        event_truth->Draw(( TString("pt >> h_pt_dir_") += i ), cut_acceptance && cut_promptphoton && arr_cut_isolation.At(i)->GetTitle() );

        cout << "Entries (All)   : " << h_pt_all->GetEntries() << endl;
        cout << "Entries (Pi0)   : " << h_pt_pi0->GetEntries() << endl;
        cout << "Entries (Prompt): " << h_pt_dir->GetEntries() << endl;

        arr_hist_pt_all.Add( h_pt_all );
        arr_hist_pt_pi0.Add( h_pt_pi0 );
        arr_hist_pt_dir.Add( h_pt_dir );
      }
  }

  /* Create legend (use same for all plots) */
  {
    TLegend* legend_c0 = new TLegend(0.7,0.7,0.89,0.89);
    legend_c0->SetFillColor(0);
    legend_c0->SetLineColor(0);
    legend_c0->AddEntry(h_pt_all,"All photons","lep");
    legend_c0->AddEntry(h_pt_pi0,"#pi^{0} photons","lep");
    legend_c0->AddEntry(h_pt_dir,"Promtp photons","lep");
  }

  /* PLOT: pT_gamma spectrum for Direct Photons and Pi0 Decay Photons and All Photons */
  {
    for ( unsigned i = 0; i < arr_cut_isolation.GetEntries(); i++ )
      {
        TCanvas *c0 = new TCanvas();
        c0->SetLogy();
        arr_hist_pt_all.At(i)->DrawClone();
        arr_hist_pt_pi0.At(i)->DrawClone("sames");
        arr_hist_pt_dir.At(i)->DrawClone("sames");
        legend_c0->Draw();
        gPad->RedrawAxis();
      }
  }

  /* Plot: ratios pass cut / total */
  {
    for ( unsigned i = 0; i < arr_cut_isolation.GetEntries(); i++ )
      {
        TH1F* h_pt_all_pass_c1 = (TH1F*) arr_hist_pt_all.At(i)->Clone("h_pt_all_pass_c1");
        TH1F* h_pt_pi0_pass_c1 = (TH1F*) arr_hist_pt_pi0.At(i)->Clone("h_pt_pi0_pass_c1");
        TH1F* h_pt_dir_pass_c1 = (TH1F*) arr_hist_pt_dir.At(i)->Clone("h_pt_dir_pass_c1");

        h_pt_all_pass_c1->Divide( (TH1*)( arr_hist_pt_all.At(0) ) );
        h_pt_pi0_pass_c1->Divide( (TH1*)( arr_hist_pt_pi0.At(0) ) );
        h_pt_dir_pass_c1->Divide( (TH1*)( arr_hist_pt_dir.At(0) ) );

        h_pt_all_pass_c1->GetYaxis()->SetTitle("fraction passed cut");
        h_pt_pi0_pass_c1->GetYaxis()->SetTitle("fraction passed cut");
        h_pt_dir_pass_c1->GetYaxis()->SetTitle("fraction passed cut");

        TCanvas *c2 = new TCanvas();
        h_pt_all_pass_c1->DrawClone("");
        h_pt_pi0_pass_c1->DrawClone("sames");
        h_pt_dir_pass_c1->DrawClone("sames");
        legend_c0->Draw();
        gPad->RedrawAxis();
      }
  }

  return 0;
}
