#include "EmcLocalRecalibratorSasha.h"

#include <PhotonContainer.h>
#include <Photon.h>
#include <PhotonERT.h>

#include <cmath>

using namespace std;

EmcLocalRecalibratorSasha::EmcLocalRecalibratorSasha()
{
  for( int i=0; i<NMAXTWR; i++ )
  {
    fCorrTof[i]=0.;
  }
}

void EmcLocalRecalibratorSasha::ApplyClusterCorrection( const int runno, PhotonContainer *photoncont )
{
  for ( unsigned i = 0; i < photoncont->Size(); i++ )
  {
    Photon *photon = photoncont->GetPhoton(i);
    float eclcore = photon->get_Ecorr();
    float emc_tof = photon->get_tof();
    int id = photon->get_towerid();

    if( eclcore > 0.1 ) {
      emc_tof -= GetTowerTofCorr(runno,id,eclcore);
    }

    photoncont->GetPhoton(i)->set_E( eclcore );
    photoncont->GetPhoton(i)->set_tof( emc_tof );
  }

  return;
}

void EmcLocalRecalibratorSasha::anaGetCorrTof(const char* fname)
{
  FILE *pf;

  // Tower-by-tower corrections

  pf = fopen(fname,"r");
  if( !pf ) {
    printf("Error in GetCorrTof: error opening tower ToF correction file %s\n",fname);
    printf("No tower ToF correction will be applied\n");
  }
  else {

    float corr;
    int ich, code;
    int nn=0;
    for( int i=0; i<NMAXTWR; i++ ) {

      if( fscanf(pf, "%d %d %f\n", &ich, &code, &corr) != 3 ) {
        printf("Error in GetCorrTof: error reading data from %s\n",fname);
        return;
      }

      // Attention!! I use my indexing for warn map in this program!!!
      if( ich >= 10368 && ich < 15552 ) { // PbSc
        if( ich < 12960 ) ich += 2592;
        else              ich -= 2592;
      }
      else if( ich >= 15552 )           { // PbGl
        if( ich < 20160 ) ich += 4608;
        else              ich -= 4608;
      }

      fCorrTof[ich] = corr;
      nn++;
    }
  
    fclose(pf);
    printf("Info from GetCorrTof: %d towers read from file %s\n",nn,fname);
  
  }
}

float EmcLocalRecalibratorSasha::GetTowerTofCorr(int run, int id, float en)
{
  int sec;
  float t1, t2;

  if( id<0 || id>24767 ) return 0;
  t1 = fCorrTof[id];

  if( id<15552 ) sec = id/2592;
  else           sec = 6 + (id-15552)/4608;

  //  int irun = FindRun_db(run);
  //  if( irun < 0 ) t2 = 0;
  //  else t2 = fCorrTofRun[irun][sec];

  // Temporary solution for run dependence in Run13
  //
  const int NSEC = 8;
  const int NPERIOD = 4;
  int run_min[NPERIOD] = {386700,391813,394266,395223};
  int run_max[NPERIOD] = {391812,394265,395222,398200};
  float tshift[NSEC][NPERIOD] = { // This is from calt_run.C
    {-8.063,   -3.896,   -3.970,   -9.230},
    {-8.350,   -4.353,   -4.419,   -9.702},
    {-13.340,   -8.254,   -8.245,   -9.165},
    {-13.254,   -8.107,   -8.098,   -9.022},
    {-9.798,   -8.608,   -8.513,   -9.034},
    {-9.811,   -8.771,   -8.850,   -9.010},
    {4.738,   -4.713,    4.016,    3.719},
    {4.992,   -4.714,    4.292,    3.985}
  };

  int iper = 0;
  for( int ip=0; ip<NPERIOD; ip++ ) {
    if( run>=run_min[ip] && run<=run_max[ip] ) {
      iper = ip;
      continue;
    }
  }
  t2 = tshift[sec][iper]-tshift[sec][0];

  // Now energy (slewing) correction in Run13 from calt.C::global()
  float dt = 0;
  if( en>0.1 && id<15552 ) { // Only for PbSc
    dt = -13./pow(en,0.45)+7.79; // Only for Run13 !!!!!
    if( dt<-15 ) dt = -15;
  }

  return t1+t2+dt;
}
