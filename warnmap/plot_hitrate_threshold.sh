plot_hitrate_threshold( string checkfile="warnmap-output/Checkplots_GernerateWarnmap_nsigma10_niter10_erange2.root" , bool writeplots = true )
{
  gStyle->SetOptStat(0);

  /* Open file */
  TFile *fcheck = new TFile( checkfile.c_str(), "OPEN" );

  /* Loop over sectors */
  for ( int sector = 0; sector < 8; sector++ )
    {
