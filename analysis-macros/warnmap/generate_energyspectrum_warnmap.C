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


int generate_energyspectrum_warnmap( string histfile="", string histname="none", string warnmapfile="warnmap-final/Warnmap_Run13pp510MinBias_Final.txt" , string basename_plots="plots/energyspectrum_", bool writeplots = true )
{
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

  /* create new histograms for each sector */
  float energybins[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
		     7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
		     13.0, 14.0, 15.0 };
  TH1D* h_base = new TH1D("h_base","",15,energybins);
  h_base->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_base->GetYaxis()->SetTitle("# entries / 1e6");
  //  h_base->GetYaxis()->SetRangeUser(?,?);

  TH1D* h_energy_sector_nowarn[8];
  h_energy_sector_nowarn[0] = (TH1D*)h_base->Clone("energy_sector_0_nowarn");
  h_energy_sector_nowarn[1] = (TH1D*)h_base->Clone("energy_sector_1_nowarn");
  h_energy_sector_nowarn[2] = (TH1D*)h_base->Clone("energy_sector_2_nowarn");
  h_energy_sector_nowarn[3] = (TH1D*)h_base->Clone("energy_sector_3_nowarn");
  h_energy_sector_nowarn[4] = (TH1D*)h_base->Clone("energy_sector_4_nowarn");
  h_energy_sector_nowarn[5] = (TH1D*)h_base->Clone("energy_sector_5_nowarn");
  h_energy_sector_nowarn[6] = (TH1D*)h_base->Clone("energy_sector_6_nowarn");
  h_energy_sector_nowarn[7] = (TH1D*)h_base->Clone("energy_sector_7_nowarn");

  h_energy_sector_nowarn[0]->SetLineColor(kOrange+5);
  h_energy_sector_nowarn[1]->SetLineColor(kGray+2);
  h_energy_sector_nowarn[2]->SetLineColor(4);
  h_energy_sector_nowarn[3]->SetLineColor(6);
  h_energy_sector_nowarn[4]->SetLineColor(8);
  h_energy_sector_nowarn[5]->SetLineColor(9);
  h_energy_sector_nowarn[6]->SetLineColor(1);
  h_energy_sector_nowarn[7]->SetLineColor(2);

  TH1D* h_energy_sector[8];
  h_energy_sector[0] = (TH1D*)h_base->Clone("energy_sector_0");
  h_energy_sector[1] = (TH1D*)h_base->Clone("energy_sector_1");
  h_energy_sector[2] = (TH1D*)h_base->Clone("energy_sector_2");
  h_energy_sector[3] = (TH1D*)h_base->Clone("energy_sector_3");
  h_energy_sector[4] = (TH1D*)h_base->Clone("energy_sector_4");
  h_energy_sector[5] = (TH1D*)h_base->Clone("energy_sector_5");
  h_energy_sector[6] = (TH1D*)h_base->Clone("energy_sector_6");
  h_energy_sector[7] = (TH1D*)h_base->Clone("energy_sector_7");

  h_energy_sector[0]->SetLineColor(kOrange+5);
  h_energy_sector[1]->SetLineColor(kGray+2);
  h_energy_sector[2]->SetLineColor(4);
  h_energy_sector[3]->SetLineColor(6);
  h_energy_sector[4]->SetLineColor(8);
  h_energy_sector[5]->SetLineColor(9);
  h_energy_sector[6]->SetLineColor(1);
  h_energy_sector[7]->SetLineColor(2);

  /* Loop over input histograms and fill energy spectrum histograms */

  /* open input file */
  TFile *fhin = new TFile("warnmap-data/WarnmapData_Run13pp510MinBias.root","OPEN");

  TH2F* h_hits = (TH2F*)fhin->Get("hitmap_energy");

  /* counter for excluded towers */
  int tower_excluded[8];
  for ( int s = 0; s< 8; s++ )
    {
      tower_excluded[s] = 0;
    }

  /* loop over all towers */
  for ( int itowerid = 0; itowerid < 24768; itowerid++ )
    {
      int isector = 0;
      int izbin = 0;
      int iybin = 0;
      int istatus = 0;

      towerid2position( itowerid, isector, iybin, izbin );

      //      cout << itowerid << " " << isector << endl;
      //      cout << isector << " += " << icounts << endl;

      /* count towers that are excluded based on warnmap */
      bool iexclude = false;
      if ( warnmap[isector][iybin][izbin] > 10 )
	{
	  iexclude = true;
	  tower_excluded[isector]++;
	}

      /* loop over all energy ranges */
      for ( int ienergybin = 0; ienergybin < 20; ienergybin++ )
	{
	  float ienergy = h_hits->GetYaxis()->GetBinCenter(ienergybin);

	  float icounts = h_hits->GetBinContent( itowerid+1 , ienergybin );

	  h_energy_sector_nowarn[isector]->Fill( ienergy, icounts );

	  /* only include towers with good status on warnmap */
	  if ( ! iexclude )
	    {
	      h_energy_sector[isector]->Fill( ienergy, icounts );
	    }
	}
    }

  /* Scale down energy spectrum by 1e6 for more manageable axis range */
  for ( int s = 0; s < 8; s++ )
    {
      h_energy_sector_nowarn[s]->Scale(1./1000000.);
      h_energy_sector[s]->Scale(1./1000000.);
    }

  /* Correct energy spectrum for number of excluded towers in each sector */
  for ( int s = 0; s < 8; s++ )
    {
      float livescalefactor = ( ntower_per_sector[s] / ( ntower_per_sector[s] - (float)tower_excluded[s] ) );
      h_energy_sector[s]->Scale(livescalefactor);
      cout << "Livescalefactor: " << livescalefactor << endl;
    }

  /* Plot energy spectrum */
  h_base->SetLineColor(0);
  h_base->GetYaxis()->SetRangeUser(0.001,50000);

  TCanvas *c1 = new TCanvas();
  c1->SetLogy();
  h_base->Draw();
  for ( int s =0; s < 8; s++ )
    h_energy_sector_nowarn[s]->Draw("same");
  gPad->RedrawAxis();

  c1->Print("plots/energySpectrum_Run13pp510MB_nowarn.eps");
  c1->Print("plots/energySpectrum_Run13pp510MB_nowarn.png");

  TCanvas *c2 = new TCanvas();
  c2->SetLogy();
  h_base->Draw();
  for ( int s =0; s < 8; s++ )
    h_energy_sector[s]->Draw("same");
  gPad->RedrawAxis();

  c2->Print("plots/energySpectrum_Run13pp510MB.eps");
  c2->Print("plots/energySpectrum_Run13pp510MB.png");

  TCanvas *c3 = new TCanvas();
  h_base->Draw();
  TLegend *leg = new TLegend(0.2,0.6,0.8,0.85);
  leg->SetNColumns(2);
  leg->AddEntry(h_energy_sector[0],"sector 0","l");
  leg->AddEntry(h_energy_sector[1],"sector 1","l");
  leg->AddEntry(h_energy_sector[2],"sector 2","l");
  leg->AddEntry(h_energy_sector[3],"sector 3","l");
  leg->AddEntry(h_energy_sector[4],"sector 4","l");
  leg->AddEntry(h_energy_sector[5],"sector 5","l");
  leg->AddEntry(h_energy_sector[6],"sector 6","l");
  leg->AddEntry(h_energy_sector[7],"sector 7","l");
  leg->Draw();
  c3->Print("plots/energySpectrum_legend.eps");
  c3->Print("plots/energySpectrum_legend.png");

  return 0;
}
