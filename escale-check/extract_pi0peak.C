int extract_pi0peak( unsigned run_index, int run_number, string histfile, string histname, bool visualize )
{
  /* Get number of events in run */
  int nevents = 0;

  /* Open 2D histogram file */
  TFile *f_in = new TFile( histfile.c_str(), "OPEN" );
  TH2F* h_inv_mass_allsector = (TH2F*) f_in->Get( histname.c_str() );
  cout << "Entries (all sectors): " << h_inv_mass_allsector->GetEntries() << endl;

  /* parameters for Gaussian fit to pi0 mass peak */
  float fitrange_min = 0.10;
  float fitrange_max = 0.17;

  /* collect sector information in vectors */
  vector< int > v_sector;
  vector< float > v_mean;
  vector< float > v_dmean;
  vector< float > v_sigma;
  vector< float > v_dsigma;
  vector< float > v_chisquare;

  /* Project sector bins into single 1-D histogram */
  TH1F* h_inv_mass_sector[8];
  for ( unsigned int s = 0; s < 8; s++ )
    {
      TString hname_s("h_inv_mass_sector");
      hname_s+=s;

      //cout << "Processing " << hname_s << endl;
      h_inv_mass_sector[s] = (TH1F*)h_inv_mass_allsector->ProjectionX( hname_s, s+1, s+1, "e" ); // sector/iterator count from 0, ROOT bins count from 1

      TF1 *f_fit = new TF1("f_fit","gaus");

      (h_inv_mass_sector[s])->Rebin(7);
      (h_inv_mass_sector[s])->Fit(f_fit, "Q", "", fitrange_min, fitrange_max);

      v_sector.push_back( s );
      v_mean.push_back( f_fit->GetParameter(1) );
      v_dmean.push_back( f_fit->GetParError(1) );
      v_sigma.push_back( f_fit->GetParameter(2) );
      v_dsigma.push_back( f_fit->GetParError(2) );
      v_chisquare.push_back( f_fit->GetChisquare() );
    }

  /* print results */
  for ( unsigned int s = 0; s < 8; s++ )
    {
      cout << "pi0peak "
	   << run_index << " "
	   << run_number << " "
	   << nevents << " "
	   << v_sector.at( s ) << " "
	   << v_mean.at( s ) << " "
	   << v_dmean.at( s ) << " "
	   << v_sigma.at( s ) << " "
	   << v_dsigma.at( s ) << " "
	   << v_chisquare.at( s )
	   << endl;
    }

  /* Visual output */
  if ( visualize )
    {
      //      gStyle->SetOptStat(0);

      /* Plot energy spectrum */
      TCanvas *c1 = new TCanvas();
      //c1->SetLogy();
      h_inv_mass_sector[0]->Draw();

      //      c1->Print("plots_escale/two_photon_invariant_mass_allsector.eps");
      //      c1->Print("plots_escale/two_photon_invariant_mass_allsector.png");
    }

  return 0;
}
