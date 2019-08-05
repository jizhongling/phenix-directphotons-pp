#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

void IsGoodTree(const char *fname, const double nevents)
{
  TFile *f = new TFile(fname);
  if(!f || f->IsZombie()) exit(1);
  TTree *t = (TTree*)f->Get("T");
  if(!t || t->GetEntries() < nevents) exit(1);
  exit(0);
}

void IsGoodHisto(const char *fname, const double nevents)
{
  TFile *f = new TFile(fname);
  if(!f || f->IsZombie()) exit(1);
  TH1 *h = (TH1*)f->Get("h_events");
  if(!h || h->GetEntries() < nevents) exit(1);
  exit(0);
}

void IsGoodFile(int argc, char *argv[])
{
  if(argc < 4)
  {
    printf("Usage: %s <tree|histo> <file.root> <nevents>\n", argv[0]);
    exit(1);
  }
  double nevents = atof(argv[3]) - 0.5;
  if(strcmp(argv[1],"tree") == 0)
    IsGoodTree(argv[2],nevents);
  else if(strcmp(argv[1],"histo") == 0)
    IsGoodHisto(argv[2],nevents);
  exit(1);
}

#ifndef __CINT__
#include "TApplication.h"

void StandaloneApplication(int argc, char *argv[]) {
  IsGoodFile(argc, argv);
}

int main(int argc, char *argv[]) {
  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplication(app.Argc(), app.Argv());
  //app.Run();
  return 0;
}

#endif
