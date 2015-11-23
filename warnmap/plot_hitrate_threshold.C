plot_hitrate_threshold( string checkfile="warnmap-output/Checkplots_GernerateWarnmap_nsigma10_niter10_erange2.root" , bool writeplots = true )
{
  gStyle->SetOptStat(0);

  /* Open file */
  TFile *fcheck = new TFile( checkfile.c_str(), "OPEN" );

  /* Loop over sectors */
  for ( int sector = 0; sector < 8; sector++ )
    {
      /* Read histogram */
      TString histname("tower_hitfrequency_sector_");
      histname += sector;

      TH1I* h_rate = (TH1I*)fcheck->Get( histname );

      /* Read graph with threshold positions */
      TString graph_threshold_name("g_iteration_threshold_sector_");
      graph_threshold_name += sector;

      TGraph* g_threshold = (TGraph*)fcheck->Get( graph_threshold_name );
      int n_thresholds = g_threshold->GetN();
      double* thresh = g_threshold->GetY();

//      TLine *lines[100];
//      for ( int i = 0; i < n_thresholds; i++ )
//	{
//	  cout << "Threshold " << thresh[i] << endl;
//	  lines[i] = new TLine(thresh[i],0,thresh[i],100000);
//	  lines[i]->SetLineColor(kRed);
//	}

      /* Draw histogram with thresholds */
      TCanvas *c1 = new TCanvas();

      h_rate->Draw();

      // call c1->Update() to make Uymax available
      c1->Update();

      TLine *line_thr = new TLine(thresh[n_thresholds-1],
				  0,
				  thresh[n_thresholds-1],
				  gPad->GetUymax() );
      line_thr->SetLineColor(kRed);
      line_thr->Draw();

      c1->SetLogy(1);

      TString filename("plots/hitrate_threshold_sector_");
      filename+=sector;
      filename+=".eps";

      TString filenamep("plots/hitrate_threshold_sector_");
      filenamep+=sector;
      filenamep+=".png";

      c1->Print( filename );
      c1->Print( filenamep );
    }

  return;
}
