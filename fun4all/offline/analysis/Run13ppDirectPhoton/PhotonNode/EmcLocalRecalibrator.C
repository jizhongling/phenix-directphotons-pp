#include "EmcLocalRecalibrator.h"

#include "AnaToolsTowerID.h"
#include "AnaToolsPhoton.h"

#include <PhotonContainer.h>
#include <Photon.h>
#include <PhotonERT.h>

#include <phool.h>

#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

EmcLocalRecalibrator::EmcLocalRecalibrator() :
  _file_tofmap(""),
  _file_energycalibration("")
{
  datatype = ERT;

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

  /* Define function for removing pT dependance on PbSc */
  _pbsc_recor_func = new TF1("pbsc_recor_func", "6.528*sqrt(x-0.7062)-4.786-1.014*x+0.007968*x**2", 0.01, 100);

}

void EmcLocalRecalibrator::SelectMB()
{
  datatype = MB;
  return;
}

void EmcLocalRecalibrator::SelectERT()
{
  datatype = ERT;
  return;
}

void EmcLocalRecalibrator::ApplyClusterCorrection( PhotonContainer *photoncont )
{
  for ( unsigned i = 0; i < photoncont->Size(); i++ )
  {
    double ecore_corr = GetCorrectedEcore( photoncont->GetPhoton(i) );
    photoncont->GetPhoton(i)->set_E( ecore_corr );

    double tof_corr = GetCorrectedTof( photoncont->GetPhoton(i) );
    photoncont->GetPhoton(i)->set_tof( tof_corr );
  }
  return;
}

double EmcLocalRecalibrator::GetCorrectedEcore( const Photon *photon )
{
  int sector = anatools::GetSector(photon);

  double ecore_raw = photon->get_E();

  /* Apply non-linearity correction separately for PbScintillator and PbGlass calorimeter */
  double ecore_lin = ecore_raw;

  if ( sector < 6 )
    ecore_lin /= _pbsc_cor_func->Eval(ecore_raw);
  else
    ecore_lin /= _pbgl_cor_func->Eval(ecore_raw);

  /* Apply run-by-run offset */
  double ecore = ecore_lin * _energycalibration[ sector ];

  return ecore;
}

double EmcLocalRecalibrator::GetCorrectedTof(const Photon *photon)
{
  double tof_raw = photon->get_tof();

  int sector = anatools::GetSector(photon);
  int towerid = photon->get_towerid();

  double tofcorrection = _tofmap[towerid];

  // remove pT depence for PbSc in ERT sample
  if ( datatype == ERT && sector < 6 )
  {
    double pT = anatools::Get_pT(photon);
    tofcorrection += _pbsc_recor_func->Eval(pT);
  }

  if(tofcorrection<-999) tofcorrection = 0;

  return tof_raw - tofcorrection;
}

void EmcLocalRecalibrator::ReadEnergyCorrection(const int& a_runnumber)
{
  if( _file_energycalibration == "" )
  {
    cout << PHWHERE << "File for Tower Energy Calibration NOT SPECIFIED." << endl;
    exit(1);
  }

  ifstream calib_fin;
  calib_fin.open(_file_energycalibration.c_str());

  bool chk = false;
  double con[8];

  /* loop over lines in file */
  string calib_line;
  while ( getline( calib_fin, calib_line ) )
  {
    if( calib_line == "" ) break;

    istringstream iss(calib_line);

    int run;

    iss >> run >> con[0] >> con[1] >> con[2] >> con[3] >> con[4] >> con[5] >> con[6] >> con[7];

    if(run==a_runnumber)
    {
      chk = true;
      break;
    }
  }

  for(int i=0;i<8;++i)
  {
    int arm = i / 4;
    int rawsector = i % 4;
    int sector = anatools::CorrectClusterSector( arm, rawsector );

    if(chk==true) _energycalibration[sector] = con[i];
  }

  calib_fin.close();

  return;
}

void EmcLocalRecalibrator::ReadTofCorrection( const int& a_fillnumber )
{
  if( _file_tofmap == "" )
  {
    cout << PHWHERE << "File for TOFMap NOT SPECIFIED." << endl;
    exit(1);
  }

  cout << "Read from file " << _file_tofmap << endl;

  int fillnumber = a_fillnumber;

  TFile* tofmapf = new TFile(_file_tofmap.c_str());
  if(!tofmapf->IsOpen())
  {
    cout << PHWHERE << "Can not find tofmap file." << endl;
    delete tofmapf;
    exit(1);
  }

  TTree* Ttof = (TTree*)tofmapf->Get("T");
  if(!Ttof)
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
  bool chk=false;
  for(int i=0;i<nentries;++i)
  {
    Ttof->GetEntry(i);

    if(fill==fillnumber)
    {
      cout << "Correct TOF corrrection entry is found." << endl;
      chk = true;
      break;
    }
  }
  if(chk==false)
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
}
