plot_warnmap2D( string warnmapfile , bool writeplots = true )
{
  gStyle->SetOptStat(0);

  /* Read warnmap as tree */
  TTree *twarn = new TTree();
  twarn->ReadFile( warnmapfile.c_str(), "sector/I:y:z:status" );

  // Declaration of leaf types
  Int_t           sector;
  Int_t           y;
  Int_t           z;
  Int_t           status;

  // List of branches
  TBranch        *b_sector;   //!
  TBranch        *b_y;   //!
  TBranch        *b_z;   //!
  TBranch        *b_status;   //!

  twarn->SetBranchAddress("sector", &sector, &b_sector);
  twarn->SetBranchAddress("y", &y, &b_y);
  twarn->SetBranchAddress("z", &z, &b_z);
  twarn->SetBranchAddress("status", &status, &b_status);



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

      h_warnmap[sector]->Fill( z , y , status );

    }

  /* build base filename for plots */
  std::size_t pos = warnmapfile.find("/");
  string filename_cut = warnmapfile.substr( pos+1 );
  (filename_cut.erase( filename_cut.length()-4 ,4 ));//.erase(0,25);


  /* plot warnmaps */
  for( int sector = 0; sector < 8; sector++ )
    {
      TCanvas *c1 = new TCanvas();
      (h_warnmap[sector])->Draw("colz");

      if ( writeplots )
	{
	  TString filename("plots-map2D/warnmap2D_");
	  filename+=filename_cut;
	  filename+="_sector_";
	  filename+=sector;
	  filename+=".eps";

	  TString filenamep("plots-map2D/warnmap2D_");
	  filenamep+=filename_cut;
	  filenamep+="_sector_";
	  filenamep+=sector;
	  filenamep+=".png";

	  c1->Print( filename );
	  c1->Print( filenamep );
	}

    }


  /* Draw combined canvas with all sectors */
  TCanvas *c2 = new TCanvas();
  c2->SetCanvasSize( 1000, 500 );
  c2->SetWindowSize( 1000, 500 );
  c2->Divide(4,2);

  for( int sector = 0; sector < 8; sector++ )
    {
      c2->cd(sector+1);
      (h_warnmap[sector])->Draw("colz");
    }

  TString filenameAll("plots-map2D/warnmap2D_");
  filenameAll+=filename_cut;
  filenameAll+="_sector_all.eps";

  TString filenameAllp("plots-map2D/warnmap2D_");
  filenameAllp+=filename_cut;
  filenameAllp+="_sector_all.eps";

  c2->Print( filenameAll );
  c2->Print( filenameAllp );



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



  return;
}
