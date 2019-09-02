#include "EmcLocalRecalibrator.h"

#include <emcClusterContainer.h>
#include <emcClusterContent.h>
#include <TH1.h>

#include <iostream>
#include <sstream>
#include <string>

using namespace std;
/////////////////////////////////////////////////////////////////

EmcLocalRecalibrator::EmcLocalRecalibrator() : _file_warnmap(""),
  _file_tofmap(""),
  _file_energycalibration("")
{

  for(int i=0;i<8;i++)
  {
    _energycalibration[i] = 0;
  }

  for(int j=0;j<25000;j++)
  {
    _tofmap[j] = 0;
  }

  /* Define functions for non-linearity correction of EMcal cluster */
  _pbsc_cor_func = new TF1("pbsc_cor_func", "0.003+(1-0.010/x)", 0.01, 100);
  _pbgl_cor_func = new TF1("pbgl_cor_func", "0.021+(1-0.020/x)", 0.01, 100);

}//EmcLocalRecalibrator::EmcLocalRecalibrator


void EmcLocalRecalibrator::ApplyClusterCorrection( emcClusterContainer* data_emccontainer )
{
  for ( unsigned i = 0; i < data_emccontainer->size(); i++ )
  {
    double ecore_corr = GetCorrectedEcore( data_emccontainer->getCluster(i) );
    data_emccontainer->getCluster(i)->set_ecore( ecore_corr );

    double tof_corr = GetCorrectedTof( data_emccontainer->getCluster(i) );
    data_emccontainer->getCluster(i)->set_tofcorr( tof_corr );
  }
  return;
}//EmcLocalRecalibrator::ApplyClusterCorrection


double EmcLocalRecalibrator::GetCorrectedEcore( const emcClusterContent *cluster )
{
  int sector = anatools::CorrectClusterSector( cluster->arm() , cluster->sector() );

  double ecore_raw = cluster->ecore();

  /* Apply non-linearity correction separately for PbScintillator and PbGlass calorimeter */
  double ecore_lin = ecore_raw;

  if ( sector < 6 )
    ecore_lin /= _pbsc_cor_func->Eval(ecore_raw);
  else
    ecore_lin /= _pbgl_cor_func->Eval(ecore_raw);

  /* Apply run-by-run offset */
  double ecore = ecore_lin * _energycalibration[ sector ];

  return ecore;
}//EmcLocalRecalibrator::GetCalibConst


double EmcLocalRecalibrator::GetCorrectedTof(const emcClusterContent *cluster)
{
  double tof_raw = cluster->tofcorr();

  int sector = anatools::CorrectClusterSector( cluster->arm() , cluster->sector() );
  int ytower = cluster->iypos();
  int ztower = cluster->izpos();

  int towerid = anatools::TowerID( sector , ytower , ztower );

  double tofcorrection = _tofmap[towerid];

  if(tofcorrection<-999) tofcorrection = 0;

  return tof_raw - tofcorrection;
}//EmcLocalRecalibrator::GetTofCorrection


void EmcLocalRecalibrator::ReadEnergyCorrection(const int& a_runnumber)
{
  if( _file_energycalibration == "" )
  {
    cout << PHWHERE << "File for Tower Energy Calibration NOT SPECIFIED." << endl;
    exit(1);
  }

  ifstream calib_fin;
  calib_fin.open(_file_energycalibration.c_str());

  Bool_t chk = kFALSE;
  double con[8];

  /* loop over lines in file */
  string calib_line;
  while ( getline( calib_fin, calib_line ) )
  {
    if( calib_line == "" ) break;

    istringstream iss(calib_line);

    int run;

    //if ( !( iss >> dummys >> dummy >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z ) )
    iss >> run >> con[0] >> con[1] >> con[2] >> con[3] >> con[4] >> con[5] >> con[6] >> con[7];

    if(run==a_runnumber)
    {
      chk = kTRUE;
      break;
    }
  }

  for(int i=0;i<8;++i)
  {
    int arm = i / 4;
    int rawsector = i % 4;
    int sector = anatools::CorrectClusterSector( arm, rawsector );

    if(chk==kTRUE) _energycalibration[sector] = con[i];
  }

  return;
}//EmcLocalRecalibrator::ReadRunCalibMap


void EmcLocalRecalibrator::ReadTofCorrection( const int& a_fillnumber )
{
  if( _file_tofmap == "" )
  {
    cout << PHWHERE << "File for TOFMap NOT SPECIFIED." << endl;
    exit(1);
  }

  cout << "Read from file " << _file_tofmap<< endl;

  int fillnumber = a_fillnumber;

  TFile* tofmapf = new TFile(_file_tofmap.c_str());
  if(!tofmapf->IsOpen())
  {
    cout << PHWHERE << "Can not find tofmap file." << endl;
    delete tofmapf;
    exit(1);
  }

  TTree* Ttof = (TTree*)tofmapf->Get("T");
  if(Ttof==NULL)
  {
    cout << PHWHERE << "Can not find T in tofmap root file." << endl;
    delete tofmapf;
    exit(1);
  }

  int fill;
  float tof[8][48][96];
  Ttof->SetBranchAddress("fillnumber", &fill);
  Ttof->SetBranchAddress("tof_correction", tof);

  int nentries = Ttof->GetEntries();
  Bool_t chk=kFALSE;
  for(int i=0;i<nentries;++i)
  {
    Ttof->GetEntry(i);

    if(fill==fillnumber)
    {
      cout << "Correct TOF corrrection entry is found." << endl;
      chk = kTRUE;
      break;
    }
  }
  if(chk==kFALSE)
  {
    cout << "Can not find correct TOF Correction entry." << endl;
    cout << "No TOF Correction." << endl;

    delete tofmapf;

    return;
    //exit(1);
  }

  for(int i=0; i<8; i++)
  {
    int arm = i / 4;
    int rawsector = i % 4;
    int sector = anatools::CorrectClusterSector( arm, rawsector );

    /* Max number of towers for PbSc */
    int ymax = 35;
    int zmax = 71;

    /* Different number of towers for PbGlass */
    if ( sector > 5 )
    {
      ymax = 47;
      zmax = 95;
    }

    for(int ytower=0; ytower<=ymax; ++ytower)
    {
      for(int ztower=0; ztower<=zmax; ++ztower)
      {
        int towerid = anatools::TowerID( sector , ytower , ztower );
        _tofmap[towerid] = tof[i][ytower][ztower];
      }
    }
  }

  delete tofmapf;

  return;
}//ReaMap::ReadTofMap
