void generate_warnmap()
{

  gSystem->Load("libdppp_warnmap.so");

  int nsigma = 5;
  int niterations = 10;

  //loop over energy ranges
  for ( int erange = 0; erange < 5; erange++ )
    {
      gROOT->Reset();

      stringstream ss_plotfile;
      ss_plotfile << "warnmap_output/Checkplots_GernerateWarnmap_nsigma" << nsigma << "_niter" << niterations << "_erange" << erange << ".root";

      direct_photon_pp::GenerateWarnmap *genwarn = new direct_photon_pp::GenerateWarnmap( nsigma, ss_plotfile.str() );

      stringstream ss_histoname;
      ss_histoname << "hitmap_erange_" << erange;

      int return_code = genwarn->FillHitsFromHistogram( "warnmap_data/WarnmapData_Run13pp510MinBias.root" , ss_histoname.str() );
      if( return_code < 0 )
	return;

      genwarn->FindHotTowers( niterations );

      genwarn->GeneratePlots();

      genwarn->FiducialCutHotTowers();
      genwarn->FiducialCutSectorEdges();

      stringstream ss_warnfile;
      ss_warnfile << "warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange" << erange << ".txt";
      cout << ss_warnfile.str() << endl;
      genwarn->WriteWarnmap( ss_warnfile.str() );
    }

  return;

}
