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
  histfile="/gpfs/mnt/gpfs02/phenix/spin3/nfeege/taxi_test/keep/TreeData-Run13pp510MinBias.root";
  histfile_nowarn="/gpfs/mnt/gpfs02/phenix/spin3/nfeege/taxi_test/keep/TreeData-Run13pp510MinBias.root";

  histfile="/gpfs/mnt/gpfs02/phenix/spin3/nfeege/taxi_test/keep/TreeData-Run13pp510ERT.root";
  histfile_nowarn="/gpfs/mnt/gpfs02/phenix/spin3/nfeege/taxi_test/keep/TreeData-Run13pp510ERT.root";

  histname="pT_1cluster";
  histname_nowarn="pT_1cluster_nowarn";

  gStyle->SetOptStat(0);

  /* Determine "live scale factor" for each sector */
  float livescalefactor[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
  calc_livescalefactors( livescalefactor, warnmapfile );

  /* Open 2D histogram files */
  TFile *f_in = new TFile( histfile.c_str(), "OPEN" );
  TH2F* h_pT_allSectors = (TH2F*) f_in->Get( histname.c_str() );
  cout << "Entries (all sectors): " << h_pT_allSectors->GetEntries() << endl;

  TFile *f_in_nowarn = new TFile( histfile_nowarn.c_str(), "OPEN" );
  TH2F* h_pT_allSectors_nowarn = (TH2F*) f_in_nowarn->Get( histname_nowarn.c_str() );
  cout << "Entries (all sectors, no warnmap): " << h_pT_allSectors_nowarn->GetEntries() << endl;

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


  /* Calculate bin errors */
  for ( int s = 0; s < 8; s++ )
    {
      h_pT_sector_nowarn[s]->Sumw2();
      h_pT_sector[s]->Sumw2();
    }


  /* Correct energy spectrum for number of excluded towers in each sector */
  for ( int s = 0; s < 8; s++ )
    {
      cout << "Sector " << s << " scaled by " << livescalefactor[s] << endl;
      h_pT_sector[s]->Scale(livescalefactor[s]);
    }

  /* Calculate mean of all sectors for each pT bin */
  TH1F* h_pT_mean = (TH1F*)h_pT_sector[0]->Clone( "pT_mean" );
  for ( unsigned sector = 1; sector < 8; sector++ )
    {
      cout << h_pT_mean->GetEntries() << endl;
      h_pT_mean->Add( h_pT_sector[sector] );
    }
  h_pT_mean->Scale(1./8.);

  /* Divide distributions by mean for each sector */
  TH1F* h_pT_ratio_sector[8];
  for ( unsigned sector = 0; sector < 8; sector++ )
    {
      h_pT_ratio_sector[sector] = (TH1F*)h_pT_sector[sector]->Clone( "temp" );
      h_pT_ratio_sector[sector]->Divide( h_pT_mean );
    }

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

  /* ratio plots */
  TCanvas *c4 = new TCanvas();
  TH1F* h_base_ratio = (TH1F*)h_base->Clone("h_base_ratio");
  h_base_ratio->GetYaxis()->SetRangeUser(0.2,1.8);
  h_base_ratio->GetYaxis()->SetTitle("# entries sector / mean(all sectors)");
  h_base_ratio->Draw();

  for ( int s =0; s < 8; s++ )
    h_pT_ratio_sector[s]->Draw("same");
  gPad->RedrawAxis();

  c4->Print("plots/pTSpectrum_ratio_Run13pp510ERT.eps");
  c4->Print("plots/pTSpectrum_ratio_Run13pp510ERT.png");

  return 0;
}


void calc_livescalefactors( float* livescalefactor, string warnmapfile )
{
  /* array with number of towers per sector */
  float ntower_per_sector[8] = { 2592,
				 2592,
				 2592,
				 2592,
				 2592,
				 2592,
				 4608,
				 4608 };

  /* read in tree to get number of excluded channels from warnmap */
  TTree *twarnmap = new TTree();
  twarnmap->ReadFile(warnmapfile.c_str(),"sector:ybin:zbin:status");

  for( unsigned sector = 0; sector < 8; sector++ )
    {
      TString cut_sector( "status > 10 && sector ==" );
      cut_sector += sector;

      float tower_excluded = twarnmap->GetEntries( cut_sector );

      livescalefactor[sector] = ( ntower_per_sector[sector] / ( ntower_per_sector[sector] - tower_excluded ) );

      cout << cut_sector << " -> excluding " << tower_excluded << " out of " << ntower_per_sector[sector] << " -> livescalefactor " << livescalefactor[sector] << endl;
    }

  return;
}
