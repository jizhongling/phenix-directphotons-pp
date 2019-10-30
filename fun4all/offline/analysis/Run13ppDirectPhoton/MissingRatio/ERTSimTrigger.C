#include "ERTSimTrigger.h"

#include <AnaToolsCluster.h>

#include <TOAD.h>
#include <emcClusterContent.h>

#include <TString.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include <iostream>

using namespace std;

ERTSimTrigger::ERTSimTrigger():
  f_ertsm(nullptr),
  t_ertsm(nullptr)
{
  rnd = new TRandom3();
  if(!rnd)
  {
    cout << "No rnd" << endl;
    exit(1);
  }

  ReadTriggerInfo();
}

ERTSimTrigger::~ERTSimTrigger()
{
  delete rnd;
  delete f_ertsm;
}

bool ERTSimTrigger::Fired(const emcClusterContent *cluster)
{
  const int npT = 21;
  const double pTbin[npT+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    30.0 };

  double pT = anatools::Get_pT(cluster);
  int ipt;
  for(ipt=0; ipt<npT; ipt++)
    if(pT < pTbin[ipt+1])
      break;
  if(ipt == npT)
    ipt--;

  int evtype = 2;
  int sector = anatools::GetSector(cluster);
  int sm = anatools::GetSM(cluster);
  int part = sector + 8*sm + 8*32*evtype;

  TSQLResult *res = t_ertsm->Query("value", Form("ipt==%d&&part==%d",ipt,part));
  TSQLRow *row = res->Next();
  if(row)
  {
    TString field0 = row->GetField(0); 
    double value = field0.Atof();
    delete row;
    delete res;

    if( TMath::Finite(value) && rnd->Rndm() < value )
      return true;
    else
      return false;
  }
  else
  {
    delete res;
    return false;
  }
}

void ERTSimTrigger::ReadTriggerInfo()
{
  TOAD *toad_loader = new TOAD("MissingRatio");
  toad_loader->SetVerbosity(0);
  string file_location = toad_loader->location("ERTEff-SM.root");
  cout << "TOAD file location: " << file_location << endl;

  f_ertsm = new TFile( file_location.c_str() );
  if(!f_ertsm || f_ertsm->IsZombie())
  {
    cout << "Cannot open " << file_location << endl;
    exit(1);
  }
  t_ertsm = (TTree*)f_ertsm->Get("t1");
  if(!t_ertsm)
  {
    cout << "No t_ertsm" << endl;
    exit(1);
  }

  delete toad_loader;
  return;
}
