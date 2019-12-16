// To compile: g++ -Wall -o root2pdf root2pdf.cc `root-config --cflags --libs`
#include <cstring>
#include <iostream>

#include <TApplication.h>
#include <TFile.h>
#include <TKey.h>
#include <TCanvas.h>

using namespace std;

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    cout << "Usage: " << argv[0] << " <name of root file>" << endl;
    return 1;
  }

  TApplication app("ROOT Application", &argc, argv);

  TFile *f = new TFile(app.Argv()[1]);
  if(f->IsZombie())
  {
    cout << "Cannot open file " << app.Argv()[1] << endl;
    delete f;
    return 1;
  }

  TIter next(f->GetListOfKeys());
  TKey *key;
  while((key = (TKey*)next()))
    if(strcmp(key->GetClassName(),"TCanvas") == 0)
    {
      TCanvas *c = (TCanvas*)key->ReadObj();
      c->Print(Form("histos/%s.pdf",c->GetName()));
    }

  delete f;
  return 0;
}
