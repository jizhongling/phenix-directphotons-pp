//  General PHENIX tools
#include "GenerateWarnmap.h"

#include <AnaToolsTowerID.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1I.h>
#include <TGraph.h>
#include <cstdlib>
#include <TF1.h>

#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;


direct_photon_pp::GenerateWarnmap::GenerateWarnmap( unsigned nsigma , string output_file_plots = "" )
{
  cout << "Constructing GenerateWarnmap object" << endl;

  /* Set initial values */
  for( unsigned sector = 0; sector < n_sector_; sector++ )
    {

      for( unsigned tower_z = 0; tower_z < n_tower_z_; tower_z++ )
	{

	  for( unsigned tower_y=0; tower_y < n_tower_y_; tower_y++ )
	    {

	      hitmap_[sector][tower_z][tower_y] = 0;
	      warnmap_[sector][tower_z][tower_y] = GOOD;

	    }//y

	}//z

    }//sector

  nsigma_ = nsigma;

  generate_plots_ = false;

  if ( output_file_plots != "" )
    generate_plots_ = true;

  fout_ = NULL;

  if ( generate_plots_ )
    {
      fout_ = new TFile( output_file_plots.c_str() , "RECREATE" );
      fout_->cd();
    }

  return;
}

int direct_photon_pp::GenerateWarnmap::FillHitsFromHistogram( string input_file , string input_histogram )
{

  cout << "Call FillHits" << endl;
  cout << "Open histogram " << input_histogram << " from file " << input_file << endl;

  /* Open input file */
  TFile *fin = new TFile( input_file.c_str() , "OPEN" );

  /* Get histogram from inpout file */
  TH1I *hin = (TH1I*)fin->Get( input_histogram.c_str() );

  /* Check if input histogram retrieved successfully */
  if ( hin == NULL )
    {
      cerr << "ERROR: histogram " << input_histogram << " not found in data file "
	   << input_file << ". EXIT." << endl;
      return 1;
    }

  /* Fill hitmap with entries from input histogram */
  int tower_id, tower_sector, tower_y, tower_z,  hit_counts;

  for( int bin = 1; bin <= ntowers_; bin++ )
    {

      /* bin number is (tower ID) + 1 */
      tower_id = bin - 1;

      /* Get location of tower for given tower ID */
      anatools::TowerLocation( tower_id, tower_sector, tower_y, tower_z );

      /* Read hits for this tower from histogram */
      hit_counts = hin->GetBinContent( bin );

      /* Fill hitmap */
      hitmap_[tower_sector][tower_z][tower_y] += hit_counts;

      //      cout << "Filling sector = " << tower_sector << " , y = " << tower_y << " , z = " << tower_z << " , hits = " << hit_counts << endl;

    }//bin

  /* close input file */
  fin->Close();
  fin->Delete();

  return 0;
}


int direct_photon_pp::GenerateWarnmap::FillHitsFrom2DHistogram( string input_file , string input_histogram , unsigned ybin_min , unsigned ybin_max )
{

  cout << "Call FillHits2D" << endl;
  cout << "Open histogram " << input_histogram << " from file " << input_file << endl;

  /* Open input file */
  TFile *fin = new TFile( input_file.c_str() , "OPEN" );

  /* Get histogram from inpout file */
  TH2I *hin = (TH2I*)fin->Get( input_histogram.c_str() );

  /* Check if input histogram retrieved successfully */
  if ( hin == NULL )
    {
      cerr << "ERROR: histogram " << input_histogram << " not found in data file "
	   << input_file << ". EXIT." << endl;
      return 1;
    }

  cout << "Use cluster energy range from "
       << hin->GetYaxis()->GetBinCenter( ybin_min ) - 0.5 * hin->GetYaxis()->GetBinWidth( ybin_min )
       << " GeV to "
       << hin->GetYaxis()->GetBinCenter( ybin_max ) + 0.5 * hin->GetYaxis()->GetBinWidth( ybin_max )
       << " GeV" << endl;

  /* Fill hitmap with entries from input histogram */
  int tower_id, tower_sector, tower_y, tower_z,  hit_counts;

  for( int xbin = 1; xbin <= ntowers_; xbin++ )
    {
      /* bin number is (tower ID) + 1 */
      tower_id = xbin - 1;

      /* Get location of tower for given tower ID */
      anatools::TowerLocation( tower_id, tower_sector, tower_y, tower_z );


      /* loop over seelcted y-bin range */
      for ( unsigned ybin = ybin_min; ybin <= ybin_max; ybin++ )
	{

	  /* Read hits for this tower from histogram */
	  hit_counts = hin->GetBinContent( xbin , ybin );

	  /* Fill hitmap */
	  hitmap_[tower_sector][tower_z][tower_y] += hit_counts;

	  //      cout << "Filling sector = " << tower_sector << " , y = " << tower_y << " , z = " << tower_z << " , hits = " << hit_counts << endl;

	}//ybin

    }//bin

  /* close input file */
  fin->Close();
  fin->Delete();

  return 0;
}


int direct_photon_pp::GenerateWarnmap::ReadUncalibratedTowers( std::string input_file )
{
  cout << "Look for uncalibrated towers in " << input_file << endl;

  /* Run over input txt file and get list of uncalibrated towers in
     'raw' sector numbering scheme, i.e. W0 = 0, E0 = 4. Need to convert
     to numbering scheme used in this analysis where W0 = 0, E0 = 7 */

  /* Stream to read table from file */
  ifstream istream;

  /* Open the file, if it won't open return an error */
  if (!istream.is_open())
    {
      istream.open( input_file.c_str() );
      if(!istream)
        {
          cerr << "ERROR Failed to open file " << input_file << endl;
          exit(1);
        }
    }

  string line;

  while ( getline( istream, line ) )
    {
      istringstream iss(line);

      int sector_raw = 0;
      int ytower_raw = 0;
      int ztower_raw = 0;
      int status_raw = 0;

      /* read string- break if error */
      if ( !( iss >> sector_raw >> ytower_raw >> ztower_raw >> status_raw ) )
	{
	  cerr << "ERROR in GenerateWarnmap.cc: Failed to read line in file " << input_file << endl;
	  exit(1);
	}

      int sector = sector_raw;
      if ( sector_raw > 3 )
	sector = 7 - (sector_raw - 4);

      int ytower = ytower_raw;
      int ztower = ztower_raw;

      warnmap_[sector][ztower][ytower] = UNCALIBRATED;
    }

  return 0;
}


/* Find dead towers which are towers below minimum hit frequency threshold */
int direct_photon_pp::GenerateWarnmap::FindDeadTowers( float threshold_alive )
{

  int sector = 0;
  int ytower = 0;
  int ztower = 0;

  /* Loop over all towers and set those above threshold to status hot */
  for( int id = 0; id < ntowers_; id++ )
    {
      anatools::TowerLocation(id, sector, ytower, ztower);

      if ( hitmap_[sector][ztower][ytower] < threshold_alive )
	{
	  warnmap_[sector][ztower][ytower] = DEAD;
	}
    }

  return 0;
}


/* Calculate hot tower threshold based on standard deviation of hit frequencies in each sector */
int direct_photon_pp::GenerateWarnmap::FindHotTowers()
{
  vector< float > hit_sum;
  vector< float > hit_sum_squared;

  vector< unsigned > entries;
  vector< unsigned > new_hot_found;
  vector< unsigned > all_hot_found;

  vector< float > mean;
  vector< float > sdev;

  vector< long long > threshold_hot;

  /* Initialize array values */
  for( unsigned int sector = 0; sector < n_sector_; sector++ ){
    mean.push_back( 0 );
    sdev.push_back( 0 );
    entries.push_back( 0 );
    hit_sum.push_back( 0 );
    hit_sum_squared.push_back( 0 );
    threshold_hot.push_back( 0 );
    new_hot_found.push_back( 0 );
    all_hot_found.push_back( 0 );
  }//sector

  int sector = 0;
  int ytower = 0;
  int ztower = 0;
  long counts = 0;

  for( int id = 0; id < ntowers_; id++ )
    {
      anatools::TowerLocation(id, sector, ytower, ztower);

      counts = hitmap_[sector][ztower][ytower];

      /* Skip towers previously marked as not-good */
      if( warnmap_[sector][ztower][ytower] != GOOD )
	    continue;

      entries[sector]++;
      hit_sum[sector] += counts;
      hit_sum_squared[sector] += pow(counts, 2.0);
    }//id

  for( unsigned int sector = 0; sector < n_sector_; sector++ )
    {
      // mean
      mean[sector] = ( (double)hit_sum[sector] ) / entries[sector];

      // standard deviation
      sdev[sector] = sqrt( ( (double)hit_sum_squared[sector] ) / ( entries[sector]-1 ) -
			   ( pow( mean[sector], 2.0 ) * entries[sector] / ( entries[sector] - 1 ) ) );

      threshold_hot[sector] = mean[sector] + nsigma_ * sdev[sector];
    }


  /* Loop over all towers and set those above threshold to status hot */
  for( int id = 0; id < ntowers_; id++ )
    {
      anatools::TowerLocation(id, sector, ytower, ztower);

      if ( hitmap_[sector][ztower][ytower] > threshold_hot[sector] )
	{
	  all_hot_found[sector]++;

	  if ( warnmap_[sector][ztower][ytower] != HOT )
	    new_hot_found[sector]++;

	  warnmap_[sector][ztower][ytower] = HOT;
	}
    }

  /* Save summary of this step in vecctors */
  unsigned iteration = v_iteration_.size() + 1;
  v_iteration_.push_back( iteration );
  v_thresholds_.push_back( threshold_hot );
  v_means_.push_back( mean );
  v_sdevs_.push_back( sdev );
  v_newhot_.push_back( new_hot_found );
  v_allhot_.push_back( all_hot_found );

  /* Print summary of this step */
//  for( unsigned int sector = 0; sector < n_sector_; sector++ )
//    cout << " Sector " << sector << ": Found " << new_hot_found[sector] << " new towers and " << all_hot_found[sector] << " total towers above threshold of " <<  threshold_hot[sector] << endl;

  return 0;
}


int direct_photon_pp::GenerateWarnmap::FindHotTowers( int niterations )
{
  /* Repeat find hot tower procedure for multiple iterations */
  for ( int iter_i = 1; iter_i <= niterations; iter_i++ )
    {
      cout << "FindHotTowers iteration: " << iter_i << endl;
      FindHotTowers( );
    }

  return 0;
}


int direct_photon_pp::GenerateWarnmap::FiducialCutHotTowers()
{

  //Loop over all PbSc towers
  for(int sector = 0; sector < 6; sector++ )
    {
      for( int y = 0; y < 36; y++ )
	{
	  for( int z = 0; z < 72; z++ )
	    {
	      // check if tower is flagged HOT
	      if( warnmap_[sector][z][y]==HOT )
		{
		  // loop over surrounding towers
		  for(int dz=-1; dz<=1; dz++){
		    for(int dy=-1; dy<=1; dy++){

		      // this neighbor does not exist
		      if ( (z + dz) < 0 || (y + dy) < 0 || (z + dz) >= 72 || (y + dy) >= 36 )
			continue;

		      if( warnmap_[sector][z+dz][y+dy] == GOOD )
			warnmap_[sector][z+dz][y+dy] = FIDUCIAL_HOT;
		    }//dy
		  }//dz

		}//if hot

	    }//loop z
	}//loop y
    }//loop sector


  //Loop over all PbGl towers
  for( int sector = 6; sector < 8; sector++ )
    {
      for( int y = 0; y < 48; y++ )
	{
	  for( int z = 0; z < 96; z++ )
	    {
	      if( warnmap_[sector][z][y] == HOT )
		{
		  for( int dz = -1; dz<=1; dz++ )
		    {
		    for( int dy = -1; dy<=1; dy++ )
		      {
			// this neighbor does not exist
			if ( (z + dz) < 0 || (y + dy) < 0 || (z + dz) >= 96 || (y + dy) >= 48 )
			  continue;

			if( warnmap_[sector][z+dz][y+dy] == GOOD )
			  warnmap_[sector][z+dz][y+dy] = FIDUCIAL_HOT;
		      }//dy
		    }//dz

		}//if hot
	    }//z
	}//y
    }//sector

  return 0;
}


int direct_photon_pp::GenerateWarnmap::FiducialCutSectorEdges()
{
  // Set sector edge of n_exclude_sector towers for hard fiducial cut
  int n_exclude_sector = 2;

  // Set arm edge of n_exclude_sector towers for fiducial cut for isolation cut
  int n_exclude_arm = 10;

  // Loop over all PbSc sectors
  for( int sector = 0; sector < 6; sector++ )
    {
      // Loop over towers in this sector and set 'verticle strips' fiducial cuts around the sectors
      // and on the arms
      for( int y = 0; y < 36; y++ )
	{
	  for( int z = 0; z < 72; z++ )
	    {
	      if( warnmap_[sector][z][y] == GOOD &&
		  ( z < n_exclude_sector || z >= ( 72 - n_exclude_sector ) ) )
		warnmap_[sector][z][y] = FIDUCIAL_SECTOR;

	      if( warnmap_[sector][z][y] == GOOD &&
		  ( z < n_exclude_arm || z >= ( 72 - n_exclude_arm ) ) )
		warnmap_[sector][z][y] = FIDUCIAL_ARM;
	    }//z
	}//y vertical strips on arms end


      // Loop over towers in this sector and set 'horizontal strips' fiducial cuts around the sectors
      // In sector 0, 3 and 4: set the 'horizonatal strips' fiducial cut on the arms
      for( int z = 0; z < 72; z++ )
	{
	  for( int y = 0; y < n_exclude_arm; y++ )
	    {
	      if( ( warnmap_[sector][z][y] == GOOD || warnmap_[sector][z][y] == FIDUCIAL_ARM ) &&
		  y < n_exclude_sector )
		warnmap_[sector][z][y] = FIDUCIAL_SECTOR;

	      if( sector == 0 ) //arm edge
		{
		  if( warnmap_[sector][z][y] == GOOD )
		    warnmap_[sector][z][y] = FIDUCIAL_ARM;
		}//sector==0
	    }//y

	  for( int y = ( 36 - n_exclude_arm ); y < 36; y++ )
	    {
	      if( ( warnmap_[sector][z][y] == GOOD || warnmap_[sector][z][y] == FIDUCIAL_ARM ) &&
		  y >= 36 - n_exclude_sector )
		warnmap_[sector][z][y] = FIDUCIAL_SECTOR;

	      if( sector == 3 || sector == 4 ) //arm edge
		{
		  if( warnmap_[sector][z][y] == GOOD )
		    warnmap_[sector][z][y] = FIDUCIAL_ARM;
		}// sector == 3||4

	    }//y

	}//z horizontal strips end

    }//sector (PbSc)

  //-----

  // Loop over all PbGl sectors
  for( int sector = 6; sector < 8; sector++ )
    {
      // Loop over towers in this sector and set 'verticle strips' fiducial cuts around the sectors
      // and on the arms
      for( int y = 0; y < 48; y++ )
	{
	  for( int z = 0; z < 96; z++ )
	    {
	      if( warnmap_[sector][z][y] == GOOD &&
		  ( z < n_exclude_sector || z >= ( 96 - n_exclude_sector ) ) )
		warnmap_[sector][z][y] = FIDUCIAL_SECTOR;

	      if( warnmap_[sector][z][y] == GOOD &&
		  ( z < n_exclude_arm || z >= ( 96 - n_exclude_arm ) ) )
		warnmap_[sector][z][y] = FIDUCIAL_ARM;
	    }//z
	}//y vertical strips on arms end


      // Loop over towers in this sector and set 'horizontal strips' fiducial cuts around the sectors
      // In sector 7: set the 'horizonatal strips' fiducial cut on the arms
      for( int z = 0; z < 96; z++ )
	{
	  for( int y = 0; y < n_exclude_arm; y++ )
	    {
	      if( ( warnmap_[sector][z][y] == GOOD || warnmap_[sector][z][y] == FIDUCIAL_ARM ) &&
		  y < n_exclude_sector )
		warnmap_[sector][z][y] = FIDUCIAL_SECTOR;

	      if( sector == 7 ) //arm edge
		{
		  if( warnmap_[sector][z][y] == GOOD )
		    warnmap_[sector][z][y] = FIDUCIAL_ARM;
		}//sector==7
	    }//y

	  for( int y = ( 48 - n_exclude_arm ); y < 48; y++ )
	    {
	      if( ( warnmap_[sector][z][y] == GOOD || warnmap_[sector][z][y] == FIDUCIAL_ARM ) &&
		  y >= 48 - n_exclude_sector )
		warnmap_[sector][z][y] = FIDUCIAL_SECTOR;
	    }//y

	}//z horizontal strips end

    }//sector (PbGl)

  return 0;
}

int direct_photon_pp::GenerateWarnmap::WriteWarnmap( string filename )
{

  ofstream warnout( filename.c_str() , ios::trunc );

  int sector, ztower, ytower;

  for( int id = 0; id < ntowers_; id++ )
    {
      anatools::TowerLocation( id, sector, ytower, ztower );
      warnout << sector << " "
	      << ytower << " "
	      << ztower << " "
	      << warnmap_[sector][ztower][ytower]
	      << endl;
    }

  warnout.close();

  return 0;
}


void direct_photon_pp::GenerateWarnmap::GeneratePlots()
{
  /* Create hitrate distribution histograms for each sector */
  TH1I* h_hitrate_distribution;

  for( int sector = 0; sector < 8; sector++ )
    {
      float mean = v_means_.at( v_iteration_.size() - 1 ).at( sector );
      float sdev = v_sdevs_.at( v_iteration_.size() - 1 ).at( sector );

      double range_highend = (double)( (int)( 1.5 * nsigma_ * sdev + mean ) );
      double range_lowend = 0;

      int nbins = (int)range_highend - (int)range_lowend;

      /* scale number of bins down to a easier to handle number */
      while ( nbins > 200 )
	{
	  nbins = nbins/10;
//      if( nbins >= 20000 )
//	nbins = nbins/40;
//      if( nbins < 20000 && nbins >= 10000 )
//	nbins = nbins/20;
//      if( nbins < 10000 && nbins >= 1000 )
//	nbins = nbins/10;
//      if( nbins < 1000 && nbins >= 200 )
//	nbins = nbins/5;
	}

      cout << "Highend: " << range_highend << " , lowend: " << range_lowend << " , bins: " << nbins << endl;

      /* Create histograms */
      TString name_hitfreq("tower_hitfrequency_sector_");
      name_hitfreq += sector;
      TString title("");
      h_hitrate_distribution = new TH1I( name_hitfreq, title, nbins, range_lowend, range_highend );
      h_hitrate_distribution->GetXaxis()->SetTitle("tower hit frequency");
      h_hitrate_distribution->GetYaxis()->SetTitle("# entries");

      for(int y = 0; y < 48; y++ )
	{
	  if( sector < 6 && y > 35 )
	    continue;

	  for( int z = 0; z < 96; z++ )
	    {
	      if( sector < 6 && z > 71 )
		continue;

	      h_hitrate_distribution->Fill( hitmap_[sector][z][y] );
	    }//z
	}//y
      // h_dstrb[erng][sect]->Fit("gaus","Q","",0,nrms*rms[erng][sect] + mean[erng][sect]);
      fout_->cd();
      h_hitrate_distribution->Write();
    }//sector


  /* Create graphs for sector properties as function of hot-tower-finding-iteration for each sector */
  for( int sector = 0; sector < 8; sector++ )
    {

      TString name_threshold("g_iteration_threshold_sector_");
      name_threshold += sector;

      TString name_newhot("g_iteration_newhot_sector_");
      name_newhot += sector;

      TString name_allhot("g_iteration_allhot_sector_");
      name_allhot += sector;

      TGraph *g_threshold = new TGraph( v_iteration_.size() );
      g_threshold->SetName(name_threshold);

      TGraph *g_newhot = new TGraph( v_iteration_.size() );
      g_newhot->SetName(name_newhot);

      TGraph *g_allhot = new TGraph( v_iteration_.size() );
      g_allhot->SetName(name_allhot);

      for ( unsigned i = 0; i < v_iteration_.size(); i++ )
	{
	  g_threshold->SetPoint( i,
				 v_iteration_.at( i ),
				 v_thresholds_.at( i ).at( sector ) );

	  g_allhot->SetPoint( i,
			      v_iteration_.at( i ),
			      v_allhot_.at( i ).at( sector ) );

	  g_newhot->SetPoint( i,
			      v_iteration_.at( i ),
			      v_newhot_.at( i ).at( sector ) );
	}

      /* setting graph axis title only works after drawing the graph first */
      g_threshold->Draw("Q");
      g_threshold->SetTitle("");
      g_threshold->GetXaxis()->SetTitle("# iteration");
      g_threshold->GetYaxis()->SetTitle("threshold");

      g_newhot->Draw("Q");
      g_newhot->SetTitle("");
      g_newhot->GetXaxis()->SetTitle("# iteration");
      g_newhot->GetYaxis()->SetTitle("# new hot towers");

      g_allhot->Draw("Q");
      g_allhot->SetTitle("");
      g_allhot->GetXaxis()->SetTitle("# iteration");
      g_allhot->GetYaxis()->SetTitle("# hot towers");

      fout_->cd();
      g_threshold->Write();
      g_newhot->Write();
      g_allhot->Write();
    }

  return;
}


void direct_photon_pp::GenerateWarnmap::Statistics()
{

  return;
}


void direct_photon_pp::GenerateWarnmap::Finish()
{
  if ( fout_ )
    fout_->Close();

  return;
}
