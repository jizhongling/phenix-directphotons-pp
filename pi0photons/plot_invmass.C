int plot_invmass( string histfile="", string histname="none", bool writeplots = true )
{

  /* Default names for debugging */
  histfile="../warnmap/warnmap-data/DirectPhotonPP_Run13pp510ERT.root";
  histname="inv_mass_2photon";

  gStyle->SetOptStat(0);

  /* Open 2D histogram file */
  TFile *f_in = new TFile( histfile.c_str(), "OPEN" );
  f_in->ls();

  TH2F* h_inv_mass_allpT = (TH2F*) f_in->Get( histname.c_str() );
  cout << "Entries (all pT bins): " << h_inv_mass_allpT->GetEntries() << endl;

  h_inv_mass_allpT->Draw("colz");

  /* Info: Print all pT bin ranges */
  for ( int b = 1; b <= h_inv_mass_allpT->GetNbinsY(); b++ )
    {

      float rmin = h_inv_mass_allpT->GetYaxis()->GetBinCenter( b ) - 0.5 * h_inv_mass_allpT->GetYaxis()->GetBinWidth( b );
      float rmax = h_inv_mass_allpT->GetYaxis()->GetBinCenter( b ) + 0.5 * h_inv_mass_allpT->GetYaxis()->GetBinWidth( b );

      cout << "pT bin " << b << " ranges from " << rmin << " GeV to " << rmax << " GeV." << endl;

    }

  /* Project sector bins into single 1-D histogram */
  TH1F* h_inv_mass = (TH1F*)h_inv_mass_allpT->ProjectionX( "invMass_allpT", 1, 8 );

  /* Plot energy spectrum */
  TCanvas *c1 = new TCanvas();
  c1->SetLogy();
  h_inv_mass->Draw();

  c1->Print("plots_escale/two_photon_invariant_mass_allpT.eps");
  c1->Print("plots_escale/two_photon_invariant_mass_allpT.png");

  return 0;
}
