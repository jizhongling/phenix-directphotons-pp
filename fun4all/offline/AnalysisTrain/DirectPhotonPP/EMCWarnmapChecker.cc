#include "EMCWarnmapChecker.h"

#include "AnaToolsTowerID.h"

#include <TOAD.h>
#include <emcClusterContent.h>

#include <TFile.h>
#include <TTree.h>

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

EMCWarnmapChecker::EMCWarnmapChecker()
{
  for(int sector=0; sector<NSEC; sector++)
    for(int iypos=0; iypos<NY; iypos++)
      for(int izpos=0; izpos<NZ; izpos++)
      {
        tower_status_nils[sector][iypos][izpos] = 0;
        tower_status_sasha[sector][iypos][izpos] = 0;
        tower_status_inseok[sector][iypos][izpos] = 0;
      }

  ReadWarnmapNils();
  ReadWarnmapSasha();
  ReadWarnmapInseok();
}

bool EMCWarnmapChecker::InFiducial(const emcClusterContent *cluster)
{
  int arm = cluster->arm();
  int rawsector = cluster->sector();
  int sector = anatools::CorrectClusterSector(arm, rawsector);
  int iypos = cluster->iypos();
  int izpos = cluster->izpos();

  if( tower_status_sasha[sector][iypos][izpos] == 0 )
    return true;
  else
    return false;
}

bool EMCWarnmapChecker::IsGoodTower(const emcClusterContent *cluster)
{
  int arm = cluster->arm();
  int rawsector = cluster->sector();
  int sector = anatools::CorrectClusterSector(arm, rawsector);
  int iypos = cluster->iypos();
  int izpos = cluster->izpos();

  if( tower_status_sasha[sector][iypos][izpos] == 0 ||
      tower_status_sasha[sector][iypos][izpos] == 30 )
    return true;
  else
    return false;
}

bool EMCWarnmapChecker::IsBadTower(const emcClusterContent *cluster)
{
  int arm = cluster->arm();
  int rawsector = cluster->sector();
  int sector = anatools::CorrectClusterSector(arm, rawsector);
  int iypos = cluster->iypos();
  int izpos = cluster->izpos();

  if( tower_status_sasha[sector][iypos][izpos] > 0 &&
      tower_status_sasha[sector][iypos][izpos] < 20 )
    return true;
  else
    return false;
}

bool EMCWarnmapChecker::InFiducial(int itower)
{
  if( itower < 0 || itower >= n_twrs ) return false;
  int sector, iypos, izpos;
  anatools::TowerLocation(itower, sector, iypos, izpos);
  if( tower_status_sasha[sector][iypos][izpos] == 0 )
    return true;
  return false;
}

bool EMCWarnmapChecker::IsGoodTower(int itower)
{
  if( itower < 0 || itower >= n_twrs ) return false;
  int sector, iypos, izpos;
  anatools::TowerLocation(itower, sector, iypos, izpos);
  if( tower_status_sasha[sector][iypos][izpos] == 0 ||
      tower_status_sasha[sector][iypos][izpos] == 30 )
    return true;
  return false;
}

bool EMCWarnmapChecker::IsBadTower(int itower)
{
  if( itower < 0 || itower >= n_twrs ) return false;
  int sector, iypos, izpos;
  anatools::TowerLocation(itower, sector, iypos, izpos);
  if( tower_status_sasha[sector][iypos][izpos] > 0 &&
      tower_status_sasha[sector][iypos][izpos] < 20 )
    return true;
  return false;
}

bool EMCWarnmapChecker::PassCut(const emcClusterContent *cluster)
{
  int arm = cluster->arm();
  int rawsector = cluster->sector();
  int armsect = 4*arm + rawsector;
  int iypos = cluster->iypos();
  int izpos = cluster->izpos();

  if( tower_status_inseok[armsect][iypos][izpos] )
    return false; 
  return true;
}

int EMCWarnmapChecker::GetStatusNils(int sector, int iypos, int izpos)
{
  if( sector < 0 || sector >= NSEC ||
      iypos < 0 || iypos >= NY ||
      izpos < 0 || izpos >= NZ )
    return 9999;

  return tower_status_nils[sector][iypos][izpos];
}

int EMCWarnmapChecker::GetStatusNils(const emcClusterContent *cluster)
{
  int arm = cluster->arm();
  int rawsector = cluster->sector();
  int sector = anatools::CorrectClusterSector(arm, rawsector);
  int iypos = cluster->iypos();
  int izpos = cluster->izpos();

  return tower_status_nils[sector][iypos][izpos]; 
}

int EMCWarnmapChecker::GetStatusSasha(const emcClusterContent *cluster)
{
  int arm = cluster->arm();
  int rawsector = cluster->sector();
  int sector = anatools::CorrectClusterSector(arm, rawsector);
  int iypos = cluster->iypos();
  int izpos = cluster->izpos();

  return tower_status_sasha[sector][iypos][izpos]; 
}

void EMCWarnmapChecker::ReadWarnmapNils()
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  int sector = 0;
  int biny = 0;
  int binz = 0;
  int status = 0;

  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  toad_loader->SetVerbosity(0);
  string file_location = toad_loader->location("Warnmap_Run13pp510.txt");
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> sector >> biny >> binz >> status )
  {
    /* count tower with bad status for PbSc and PbGl */
    if ( status > 10 )
    {
      if( sector < 6 ) nBadSc++;
      else nBadGl++;
    }
    tower_status_nils[sector][biny][binz] = status;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

void EMCWarnmapChecker::ReadWarnmapSasha()
{
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  int ich = 0;
  int sector = 0;
  int biny = 0;
  int binz = 0;
  int status = 0;

  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  toad_loader->SetVerbosity(0);
  string file_location = toad_loader->location("warn_all_run13pp500gev.dat");
  cout << "TOAD file location: " << file_location << endl;
  ifstream fin( file_location.c_str() );

  while( fin >> ich >> status )
  {
    /* Attention!! I use my indexing for warn map in this program!!! */
    if( ich >= 10368 && ich < 15552 ) { // PbSc
      if( ich < 12960 ) ich += 2592;
      else              ich -= 2592;
    }
    else if( ich >= 15552 )           { // PbGl
      if( ich < 20160 ) ich += 4608;
      else              ich -= 4608;
    }

    /* Get tower location */
    anatools::TowerLocation(ich, sector, biny, binz);

    /* Count tower with bad status for PbSc and PbGl */
    if ( status > 0 )
    {
      if( sector < 6 ) nBadSc++;
      else nBadGl++;
    }
    tower_status_sasha[sector][biny][binz] = status;

    /* Mark edge towers */
    if( anatools::Edge_cg(sector, biny, binz) &&
        tower_status_sasha[sector][biny][binz] == 0 )
      tower_status_sasha[sector][biny][binz] = 20;
    /* Mark fiducial arm */
    if( anatools::ArmEdge_cg(sector, biny, binz) &&
        tower_status_sasha[sector][biny][binz] == 0 )
      tower_status_sasha[sector][biny][binz] = 30;
  }

  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl << endl;
  fin.close();
  delete toad_loader;

  return;
}

void EMCWarnmapChecker::ReadWarnmapInseok()
{
  TOAD *toad_loader = new TOAD("DirectPhotonPP");
  toad_loader->SetVerbosity(0);
  string file_location = toad_loader->location("Run13pp510_WarnMap_05.root");
  cout << "TOAD file location: " << file_location << endl;

  TFile *warnmapf = new TFile(file_location.c_str());
  if(!warnmapf->IsOpen())
  {
    cout << "Can not find " << file_location << endl;
    delete toad_loader;
    delete warnmapf;
    exit(1);
  }

  TTree *Twarn = (TTree*)warnmapf->Get("T");
  if(!Twarn)
  {
    cout << "Can not find T in warnmap root file." << endl;
    delete toad_loader;
    delete warnmapf;
    exit(1);
  }

  int warn[NSEC][NY][NZ];
  Twarn->SetBranchAddress("warnmap", warn);
  Twarn->GetEntry(0);

  for (int i=0; i<NSEC; i++)
    for (int j=0; j<NY; j++)
      for (int k=0; k<NZ; k++)
        tower_status_inseok[i][j][k] = warn[i][j][k];

  delete toad_loader;
  delete warnmapf;
  return;
}
