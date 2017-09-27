struct dataset{
  string file;
  string histogram;
};

void generate_deadmap()
{

  gSystem->Load("libdgpp_warnmap.so");

  float threshold_alive = 1;

  string histlist = "generate_warnmap_histlist.txt";

  /* Stream to read table from file */
  ifstream istream_hist;

  /* Open the file, if it won't open return an error */
  if (!istream_hist.is_open())
    {
      istream_hist.open( histlist.c_str() );
      if(!istream_hist)
	{
	  cerr << "ERROR Failed to open file " << histlist << endl;
	  exit(1);
	}
    }

  string line_hist;

  while ( getline( istream_hist, line_hist ) )
    {
      /* Skip lines starting with / including a '#' */
      if ( line_hist.find("#") != string::npos )
	{
	  cout << "SKIPPING line in file: " << line_hist << endl;
	  continue;
	}

      istringstream iss(line_hist);

      string filename = "";
      string histname = "";
      unsigned ybin_min = 0;
      unsigned ybin_max = 0;

      /* read string- break if error */
      if ( !( iss >> filename >> histname >> ybin_min >> ybin_max ) )
	    {
	      cerr << "ERROR: Failed to read line in file " << histlist << endl;
	      exit(1);
	    }

      /* Build basename for histogram and dataset */
      string filename_cut = filename;
      (filename_cut.erase( filename_cut.length()-5 ,5 )).erase(0,25);

      string histname_cut = histname;
      histname_cut.erase(0,6);

      stringstream basename;
      basename << filename_cut;
      //      basename << histname_cut;

      cout << "Processing: " << filename << " & " << histname << " & " << basename.str()
	   << " & " << ybin_min << " & " << ybin_max << endl;

      /* Reset ROOT */
      gROOT->Reset();

      /* Run code to generate warnmap */
      direct_photon_pp::GenerateWarnmap *genwarn = new direct_photon_pp::GenerateWarnmap( 0 , "deadplots.root" );

      int return_code = genwarn->FillHitsFrom2DHistogram( filename , histname.c_str() , ybin_min , ybin_max );
      if( return_code < 0 )
	return;


      genwarn->FindDeadTowers( threshold_alive );

      genwarn->ReadUncalibratedTowers( "Run13_Uncalib_Jeff.txt" );


      /* Write warnmap to text file */
      stringstream ss_warnfile;
      ss_warnfile << "warnmap-output/Deadmap_"
		  << basename.str()
		  << "_ybins" << ybin_min << "to" << ybin_max
		  << ".txt";

      cout << ss_warnfile.str() << endl;
      genwarn->WriteWarnmap( ss_warnfile.str() );

      delete genwarn;
    }
  return;
}
