#include "EmcLocalRecalibratorSasha.h"

#include <emcClusterContainer.h>
#include <emcClusterContent.h>

#include <TF1.h>

#include <fstream>
#include <string>

using namespace std;

EmcLocalRecalibratorSasha::EmcLocalRecalibratorSasha()
{
  // Init calib
  fnRun_cal = 0;

  for( int i=0; i<NRUN; i++ ) {
    runlist_cal[i] = 0;
    for( int j=0; j<8; j++ )
      fCorrCalRun[i][j] = 0;
  }

  for( int i=0; i<NMAXTWR; i++ ) {
    fCorrCal[i]=1.;
    fCorrTof[i]=0.;
  }
}


void EmcLocalRecalibratorSasha::ApplyClusterCorrection( const int runno, emcClusterContainer* data_emccontainer )
{
  float corr_run[8];
  for( int is=0; is<8; is++ ) corr_run[is] = 1.;
  int irun = FindRun(runno,runlist_cal,fnRun_cal);
  if( irun>=0 ) {
    for( int is=0; is<8; is++ ) corr_run[is] = fCorrCalRun[irun][is];
  }

  for ( unsigned i = 0; i < data_emccontainer->size(); i++ )
  {
    emcClusterContent *emccl = data_emccontainer->getCluster(i);
    float eclcore = emccl->ecore();
    float emc_tof = emccl->tofcorr();

    int id;
    int arm = emccl->arm();
    int sec = emccl->sector();
    if( arm == 1 ) sec = 7-sec;
    int iz = emccl->izpos();
    int iy = emccl->iypos();
    if( sec<6 ) id = iz + 72*iy + sec*2592;
    else id = iz + 96*iy + (sec-6)*4608 + 15552;

    if( eclcore > 0.1 ) {
      float corrcal = GetTowerCorrCal(emccl);
      corrcal *= corr_run[sec]; // Sector-by-sector correction
      eclcore *= corrcal;
      emc_tof -= GetTowerTofCorr(runno,id,eclcore);
    }

    data_emccontainer->getCluster(i)->set_ecore( eclcore );
    data_emccontainer->getCluster(i)->set_tofcorr( emc_tof );
  }

  return;
}

void EmcLocalRecalibratorSasha::anaGetCorrCal(const char* fname)
{
  FILE *pf;
  pf = fopen(fname,"r");
  if( !pf ) {
    printf("Error in GetCorrCal: error opening file %s\n",fname);
    return;
  }

  int ich, code;
  float corr, dumm;
  int nn=0;
  for( int i=0; i<NMAXTWR; i++ ) {
    if( fscanf(pf, "%d %d %f %f\n", &ich, &code, &dumm, &corr) != 4 ) {
      printf("Error in GetCorrCal: error reading data from %s\n",fname);
      fclose(pf);
      return;
    }

    //    if( ich != i ) {
    //      printf("Error in GetCorrCal: wrong data format %d != %d\n",ich,i);
    //     return;
    //   }

    // Attention!! I use my indexing!!!
    if( ich >= 10368 && ich < 15552 ) { // PbSc
      if( ich < 12960 ) ich += 2592;
      else              ich -= 2592;
    }
    else if( ich >= 15552 )           { // PbGl
      if( ich < 20160 ) ich += 4608;
      else              ich -= 4608;
    }
    fCorrCal[ich] = 1.;
    if( code==0 ) fCorrCal[ich] = corr;
    nn++;

  }

  fclose(pf);
  printf("Info from GetCorrCal: %d towers read from file %s\n",nn,fname);

  return;
}

void EmcLocalRecalibratorSasha::anaGetCorrCal_run(const char* fname)
{
  fnRun_cal = 0;

  ifstream inFile(fname);
  if(!inFile){
    printf("Error in GetCorrCal_fill: error opening file %s\n",fname);
    return;
  }

  char stmp[256];
  int irun;
  float corr0, corr1, corr2, corr3, corr4, corr5, corr6, corr7;

  while(inFile){

    inFile.getline(stmp, 255);
    if(stmp[0]=='#' || stmp[0] == 0) continue; // strip comment or empty lines
    //    printf("%s\n",str);
    sscanf(stmp,"%6d %f %f %f %f %f %f %f %f",&irun,&corr0,&corr1,&corr2,&corr3,&corr4,&corr5,&corr6,&corr7);
    runlist_cal[fnRun_cal] = irun;
    fCorrCalRun[fnRun_cal][0] = corr0;
    fCorrCalRun[fnRun_cal][1] = corr1;
    fCorrCalRun[fnRun_cal][2] = corr2;
    fCorrCalRun[fnRun_cal][3] = corr3;
    fCorrCalRun[fnRun_cal][4] = corr4;
    fCorrCalRun[fnRun_cal][5] = corr5;
    fCorrCalRun[fnRun_cal][6] = corr6;
    fCorrCalRun[fnRun_cal][7] = corr7;
    //    printf("%d %f %f %f\n",irun,corr0,corr1,corr2);
    fnRun_cal++;

  }

  inFile.close();

  printf("Info from GetCorrCal_run: %d runs read from file %s\n",fnRun_cal,fname);
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

int EmcLocalRecalibratorSasha::FindRun(int run, int* rlist, int nn)
{
  if( nn <= 0 ) return -1;

  for( int i=0; i<nn; i++ ) {
    if( run == rlist[i] ) return i;
  }
  return -1;
}

float EmcLocalRecalibratorSasha::GetTowerCorrCal(const emcClusterContent* emccl)
{
  int nhit;
  int ich;
  float en, en0;
  float etwr;

  //  static TF1 fthresh_sc("fthresh_sc","0.003+(1-0.010/x)",0.01,100);
  //  static TF1 fthresh_gl("fthresh_gl","0.021+(1-0.020/x)",0.01,100);
  //  static TF1 fthresh_sc("fthresh_sc","0.003+(1-0.020/x)",0.01,100);
  //  static TF1 fthresh_gl("fthresh_gl","0.021+(1-0.030/x)",0.01,100);
  static TF1 fthresh_sc("fthresh_sc","(1-0.020/x)",0.01,100); // for 500 GeV
  static TF1 fthresh_gl("fthresh_gl","(1-0.030/x)",0.01,100); // for 500 GeV

  static float cgl_sc = 1.; // Global energy scale correction
  static float cgl_gl = 1.; // Global energy scale correction

  float eold = 0;
  float enew = 0;
  en0 = 0;
  nhit = emccl->twrhit();
  //  if( nhit > 4 ) nhit=4; // Only 4 max towers - almost "ecore" definition
  for( int itw=0; itw<nhit; itw++ ) {
    ich = emccl->towerid(itw);
    // Attention!! I use my indexing!!!
    if( ich >= 10368 && ich < 15552 ) { // PbSc
      if( ich < 12960 ) ich += 2592;
      else              ich -= 2592;
    }
    else if( ich >= 15552 )           { // PbGl
      if( ich < 20160 ) ich += 4608;
      else              ich -= 4608;
    }
    en = emccl->partesum(itw);
    etwr = en - en0;
    en0 = en;
    eold += etwr;
    enew += (etwr*fCorrCal[ich]); // Calib. corr.
  }
  // Nonlinearity due to threshold correction
  float etold = emccl->ecore();
  float etnew = enew/eold*etold;
  if( emccl->towerid(0) < 15552 ) {
    etnew *= cgl_sc;
    etnew /= fthresh_sc.Eval(etnew);
  }
  else {
    etnew *= cgl_gl;
    etnew /= fthresh_gl.Eval(etnew);
  }

  return etnew/etold;
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
