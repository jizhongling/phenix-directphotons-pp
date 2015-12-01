void towerid2position( int towerid , int &sector , int &y , int &z )
{
  int itwr=0;

  if(towerid<0 || towerid>24767){
    cout<<"ABORT: Bad tower towerid "<<towerid<<endl;
    armsect = -1;
    ztwr = -1;
    ytwr = -1;
    return;
  }

  if(towerid<15552)
    { // pbsc
      sector=towerid/2592;
      itwr=towerid%2592;
      z=itwr%72;
      y=itwr/72;
    }
  else
    { // pbgl
      sector=6+(towerid-15552)/4608;
      itwr=(towerid-15552)%4608;
      z=itwr%96;
      y=itwr/96;
    }

  return;
}


int plot_warnmap2D_paul( string warnmapfile="warnmap-paul/iter10_rms10/warn_Run9pp500MinBias.txt" , string basename_plots="plots-map2D/paul_Run9pp500MB_warnmap2D_sector_", bool writeplots = true )
{
  gStyle->SetOptStat(0);

  /* Read warnmap as tree */
  TTree *twarn = new TTree();
  twarn->ReadFile( warnmapfile.c_str(), "towerid/I:status" );

  // Declaration of leaf types
  Int_t           towerid;
  Int_t           status;

  // List of branches
  TBranch        *b_towerid;   //!
  TBranch        *b_status;   //!

  twarn->SetBranchAddress("towerid", &towerid, &b_towerid);
  twarn->SetBranchAddress("status", &status, &b_status);

  Int_t sector = 0;
  Int_t z = 0;
  Int_t y = 0;


  /* Create 2D histograms to visualize warnmap */
  TH2I* h_warnmap[8];
  int zbins = 0;
  int ybins = 0;

  /* Create histogram for each sector */
  for( int sector = 0; sector < 8; sector++ )
    {
      if(sector <6) // PbSc sector
	{
	  zbins = 72;
	  ybins = 36;
	}
      else // PbGl sector
	{
	  zbins = 96;
	  ybins = 48;
	}

      TString name("h_warnmap_sector_");
      name += sector;

      TString title("Warnmap sector ");
      title += sector;

      cout << name <<endl;
      cout << title <<endl;

      h_warnmap[sector] = new TH2I( name, title, zbins, 0, zbins, ybins, 0, ybins );
      h_warnmap[sector]->GetXaxis()->SetTitle("z [cm]");
      h_warnmap[sector]->GetYaxis()->SetTitle("y [cm]");

    } //sector

  /* Loop over warnmap and fill status of towers into map */
  Long64_t nentries = twarn->GetEntriesFast();

  for ( Long64_t jentry = 0; jentry < nentries; jentry++ )
    {
      twarn->GetEntry(jentry);

      //cout << sector << " " << y << " " << z << " " << status << endl;
      towerid2position( towerid, sector, y, z );

      h_warnmap[sector]->Fill( z , y , status );

    }


  for( int sector = 0; sector < 8; sector++ )
    {
      TCanvas *c1 = new TCanvas();
      (h_warnmap[sector])->Draw("colz");

      if ( writeplots )
	{
	  TString filename(basename_plots.c_str());
	  filename+=sector;
	  filename+=".eps";

	  TString filenamep(basename_plots.c_str());
	  filenamep+=sector;
	  filenamep+=".png";

	  c1->Print( filename );
	  c1->Print( filenamep );
	}

    }

  //  h_warnmap_[sect]->Fill(z, y, warnmap_[sect][z][y]);


//    for(int y=0; y<48; y++)
//{
//      if(sect<6 && y>35) continue;
//      for(int z=0; z<96; z++)
//{
//	if(sect<6 && z>71) continue;
//
//	id = anatools::TowerID(sect, y, z);
//	for(int erng=0; erng<6; erng++)
//{
//	  if(hotmap_[erng][sect][z][y] == 100) continue;
//	  h_counts_hot[erng]->Fill(id, hitmap_[erng][sect][z][y]);
//	}//erng
//      }//z
//    }//y
//    h_warnmap_[sect]->Write();
//  }//sect



  return 0;
}
