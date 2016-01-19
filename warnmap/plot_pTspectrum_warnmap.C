void towerid2position( int towerid , int &sector , int &y , int &z )
{
  int itwr=0;

  if(towerid<0 || towerid>24767){
    cout<<"ABORT: Bad tower towerid "<<towerid<<endl;
    sector = -1;
    y = -1;
    z = -1;
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


int plot_pTspectrum_warnmap( string histfile="", string histname="none", string histfile_nowarn="", string histname_nowarn="none", string warnmapfile="warnmap-final/Warnmap_Run13pp510MinBias_Final.txt", bool writeplots = true )
{

  /* Default names for debugging */
  histfile="/gpfs/mnt/gpfs02/phenix/spin3/nfeege/taxi_test/keep/TreeData_Warnmap_pT.root";
  histfile_nowarn="/gpfs/mnt/gpfs02/phenix/spin3/nfeege/taxi_test/keep/TreeData_Warnmap_pT.root";

  histname="pT_1cluster";
  histname_nowarn="pT_1cluster_nowarn";


  gStyle->SetOptStat(0);

  /* array with number of towers per sector */
  float ntower_per_sector[8] = { 2592,
				 2592,
				 2592,
				 2592,
				 2592,
				 2592,
				 4608,
				 4608 };

  /* Set initial values for warnmap[sector][ybin][zbin]*/
  long warnmap[8][48][96];

  for( unsigned sector = 0; sector < 8; sector++ )
    {
      for( unsigned tower_y=0; tower_y < 48; tower_y++ )
	{
	  for( unsigned tower_z = 0; tower_z < 96; tower_z++ )
	    {
	      warnmap[sector][tower_y][tower_z] = 0;
            }//z
	}//y
    }//sector

  /* Read warnmap into array */

  /* Stream to read table from file */
  ifstream istream_mapping;

  /* Open the datafile, if it won't open return an error */
  if (!istream_mapping.is_open())
    {
      istream_mapping.open( warnmapfile.c_str() );
      if(!istream_mapping)
        {
          cerr << "ERROR Failed to open warnmap file " <<	\
	    warnmapfile << endl;
          exit(1);
        }
    }

  string line_map;

  while ( getline( istream_mapping, line_map ) )
    {
      istringstream iss(line_map);

      unsigned sector, ybin, zbin, status;

      if ( !( iss >> sector >> ybin >> zbin >> status ) )
	{
	  cerr << "ERROR Failed to read line in mapping file" << endl;
	  exit(1);
	}

      warnmap[sector][ybin][zbin] = status;
    }

  cout << warnmap[7][10][15] << endl;

  /* Open 2D histogram files */
  TFile *f_in = new TFile( histfile.c_str(), "OPEN" );
  TH2F* h_pT_allSectors = (TH2F*) f_in->Get( histname.c_str() );
  cout << "Entries (all sectors): " << h_pT_allSectors->GetEntries() << endl;

  TFile *f_in_nowarn = new TFile( histfile_nowarn.c_str(), "OPEN" );
  TH2F* h_pT_allSectors_nowarn = (TH2F*) f_in_nowarn->Get( histname_nowarn.c_str() );
  cout << "Entries (all sectors, no warnmap): " << h_pT_allSectors_nowarn->GetEntries() << endl;

  /* create individual histogram for each sector */
//  float energybins[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
//		     7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
//		     13.0, 14.0, 15.0 };
//  TH1F* h_base = new TH1D("h_base","",15,energybins);
  TH1F* h_base = (TH1F*)h_pT_allSectors_nowarn->ProjectionX( "h_base", 1, 1 );
  h_base->Clear();
  h_base->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_base->GetYaxis()->SetTitle("# entries");
  h_base->GetXaxis()->SetRangeUser(3,14);
  //  h_base->GetYaxis()->SetRangeUser(1e-4,1.0);
  h_base->GetYaxis()->SetRangeUser(1e1,5e5);
  h_base->SetLineColor(0);

  TH1F* h_pT_sector_nowarn[8];
  h_pT_sector_nowarn[0] = (TH1F*)h_pT_allSectors_nowarn->ProjectionX( "pT_sector0_nowarn", 1, 1 );
  h_pT_sector_nowarn[1] = (TH1F*)h_pT_allSectors_nowarn->ProjectionX( "pT_sector1_nowarn", 2, 2 );
  h_pT_sector_nowarn[2] = (TH1F*)h_pT_allSectors_nowarn->ProjectionX( "pT_sector2_nowarn", 3, 3 );
  h_pT_sector_nowarn[3] = (TH1F*)h_pT_allSectors_nowarn->ProjectionX( "pT_sector3_nowarn", 4, 4 );
  h_pT_sector_nowarn[4] = (TH1F*)h_pT_allSectors_nowarn->ProjectionX( "pT_sector4_nowarn", 5, 5 );
  h_pT_sector_nowarn[5] = (TH1F*)h_pT_allSectors_nowarn->ProjectionX( "pT_sector5_nowarn", 6, 6 );
  h_pT_sector_nowarn[6] = (TH1F*)h_pT_allSectors_nowarn->ProjectionX( "pT_sector6_nowarn", 7, 7 );
  h_pT_sector_nowarn[7] = (TH1F*)h_pT_allSectors_nowarn->ProjectionX( "pT_sector7_nowarn", 8, 8 );

  h_pT_sector_nowarn[0]->SetLineColor(kOrange+5);
  h_pT_sector_nowarn[1]->SetLineColor(kGray+2);
  h_pT_sector_nowarn[2]->SetLineColor(4);
  h_pT_sector_nowarn[3]->SetLineColor(6);
  h_pT_sector_nowarn[4]->SetLineColor(8);
  h_pT_sector_nowarn[5]->SetLineColor(9);
  h_pT_sector_nowarn[6]->SetLineColor(1);
  h_pT_sector_nowarn[7]->SetLineColor(2);


  TH1F* h_pT_sector[8];
  h_pT_sector[0] = (TH1F*)h_pT_allSectors->ProjectionX( "pT_sector0", 1, 1 );
  h_pT_sector[1] = (TH1F*)h_pT_allSectors->ProjectionX( "pT_sector1", 2, 2 );
  h_pT_sector[2] = (TH1F*)h_pT_allSectors->ProjectionX( "pT_sector2", 3, 3 );
  h_pT_sector[3] = (TH1F*)h_pT_allSectors->ProjectionX( "pT_sector3", 4, 4 );
  h_pT_sector[4] = (TH1F*)h_pT_allSectors->ProjectionX( "pT_sector4", 5, 5 );
  h_pT_sector[5] = (TH1F*)h_pT_allSectors->ProjectionX( "pT_sector5", 6, 6 );
  h_pT_sector[6] = (TH1F*)h_pT_allSectors->ProjectionX( "pT_sector6", 7, 7 );
  h_pT_sector[7] = (TH1F*)h_pT_allSectors->ProjectionX( "pT_sector7", 8, 8 );


  h_pT_sector[0]->SetLineColor(kOrange+5);
  h_pT_sector[1]->SetLineColor(kGray+2);
  h_pT_sector[2]->SetLineColor(4);
  h_pT_sector[3]->SetLineColor(6);
  h_pT_sector[4]->SetLineColor(8);
  h_pT_sector[5]->SetLineColor(9);
  h_pT_sector[6]->SetLineColor(1);
  h_pT_sector[7]->SetLineColor(2);

//  /* Loop over input histograms and fill energy spectrum histograms */
//
//  /* open input file */
//  TFile *fhin = new TFile("warnmap-data/WarnmapData_Run13pp510MinBias.root","OPEN");
//
//  TH2F* h_hits = (TH2F*)fhin->Get("hitmap_energy");
//
//  /* counter for excluded towers */
//  int tower_excluded[8];
//  for ( int s = 0; s< 8; s++ )
//    {
//      tower_excluded[s] = 0;
//    }
//
//  /* loop over all towers */
//  for ( int itowerid = 0; itowerid < 24768; itowerid++ )
//    {
//      int isector = 0;
//      int izbin = 0;
//      int iybin = 0;
//      int istatus = 0;
//
//      towerid2position( itowerid, isector, iybin, izbin );
//
//      //      cout << itowerid << " " << isector << endl;
//      //      cout << isector << " += " << icounts << endl;
//
//      /* count towers that are excluded based on warnmap */
//      bool iexclude = false;
//      if ( warnmap[isector][iybin][izbin] > 10 )
//	{
//	  iexclude = true;
//	  tower_excluded[isector]++;
//	}
//
//      /* loop over all energy ranges */
//      for ( int ienergybin = 0; ienergybin < 20; ienergybin++ )
//	{
//	  float ienergy = h_hits->GetYaxis()->GetBinCenter(ienergybin);
//
//	  float icounts = h_hits->GetBinContent( itowerid+1 , ienergybin );
//
//	  h_energy_sector_nowarn[isector]->Fill( ienergy, icounts );
//
//	  /* only include towers with good status on warnmap */
//	  if ( ! iexclude )
//	    {
//	      h_energy_sector[isector]->Fill( ienergy, icounts );
//	    }
//	}
//    }
//

//  /* Scale down energy spectrum by 1e6 for more manageable axis range */
//  for ( int s = 0; s < 8; s++ )
//    {
//      h_pT_sector_nowarn[s]->Scale(1./1000000.);
//      h_pT_sector[s]->Scale(1./1000000.);
//    }

  /* Calculate bin errors */
  for ( int s = 0; s < 8; s++ )
    {
      h_pT_sector_nowarn[s]->Sumw2();
      h_pT_sector[s]->Sumw2();
    }

//  /* Normalize pT spectra */
//  for ( int s = 0; s < 8; s++ )
//    {
//      h_pT_sector_nowarn[s]->Scale(1./ (h_pT_sector_nowarn[s])->Integral(2,19));
//      h_pT_sector[s]->Scale(1./ (h_pT_sector[s])->Integral(2,19) );
//    }

//  /* Correct energy spectrum for number of excluded towers in each sector */
//  for ( int s = 0; s < 8; s++ )
//    {
//      float livescalefactor = ( ntower_per_sector[s] / ( ntower_per_sector[s] - (float)tower_excluded[s] ) );
//      h_energy_sector[s]->Scale(livescalefactor);
//      cout << "Livescalefactor: " << livescalefactor << endl;
//    }
//
  /* Plot energy spectrum */
  TCanvas *c1 = new TCanvas();
  c1->SetLogy();
  h_base->Draw();
  for ( int s =0; s < 8; s++ )
    h_pT_sector_nowarn[s]->Draw("same");
  gPad->RedrawAxis();

  c1->Print("plots/pTSpectrum_Run13pp510ERT_nowarn.eps");
  c1->Print("plots/pTSpectrum_Run13pp510ERT_nowarn.png");

  TCanvas *c2 = new TCanvas();
  c2->SetLogy();
  h_base->Draw();
  for ( int s =0; s < 8; s++ )
    h_pT_sector[s]->Draw("same");
  gPad->RedrawAxis();

  c2->Print("plots/pTSpectrum_Run13pp510ERT.eps");
  c2->Print("plots/pTSpectrum_Run13pp510ERT.png");

  TCanvas *c3 = new TCanvas();
  h_base->Draw();
  TLegend *leg = new TLegend(0.2,0.6,0.8,0.85);
  leg->SetNColumns(2);
  leg->AddEntry(h_pT_sector[0],"sector 0","l");
  leg->AddEntry(h_pT_sector[1],"sector 1","l");
  leg->AddEntry(h_pT_sector[2],"sector 2","l");
  leg->AddEntry(h_pT_sector[3],"sector 3","l");
  leg->AddEntry(h_pT_sector[4],"sector 4","l");
  leg->AddEntry(h_pT_sector[5],"sector 5","l");
  leg->AddEntry(h_pT_sector[6],"sector 6","l");
  leg->AddEntry(h_pT_sector[7],"sector 7","l");
  leg->Draw();
  c3->Print("plots/pTSpectrum_legend.eps");
  c3->Print("plots/pTSpectrum_legend.png");

  return 0;
}
