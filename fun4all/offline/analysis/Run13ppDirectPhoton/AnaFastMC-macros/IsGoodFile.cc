// To compile: g++ -Wall -o IsGoodFile IsGoodFile.cc `root-config --cflags --libs`
#include <iostream>
#include <sstream>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

using namespace std;

void IsGoodTree(const char *fname, const double nevents)
{
  TFile *f = new TFile(fname);
  if(!f || f->IsZombie()) exit(1);
  TTree *t = (TTree*)f->Get("T");
  if(!t || abs(t->GetEntries()-nevents) > 0.5) exit(1);
  exit(0);
}

void IsGoodHisto(const char *fname, const double nevents)
{
  TFile *f = new TFile(fname);
  if(!f || f->IsZombie()) exit(1);
  TH1 *h = (TH1*)f->Get("h_events");
  if(!h || abs(h->GetEntries()-nevents) > 0.5) exit(1);
  exit(0);
}

int main(int argc, char *argv[])
{
  if(argc < 4)
  {
    cout << "Usage: " << argv[0] << " <tree|histo> <file.root> <nevents>" << endl;
    exit(1);
  }

  string type;
  string fname;
  double nevents;
  stringstream ss;
  ss << argv[1] << ' ' << argv[2] << ' ' << argv[3];
  ss >> type >> fname >> nevents;

  if(type.compare("tree") == 0)
    IsGoodTree(fname.c_str(), nevents);
  if(type.compare("histo") == 0)
    IsGoodHisto(fname.c_str(), nevents);
  exit(1);

  return 0;
}
