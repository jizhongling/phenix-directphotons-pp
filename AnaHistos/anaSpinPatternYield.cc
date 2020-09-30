// To compile: g++ -Wall -o anaSpinPatternYield anaSpinPatternYield.cc `root-config --cflags --libs`
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

using namespace std;

int main()
{
  /* pT bins for ALL */
  const int npT_pol = 15;

  TFile *f_rlum = new TFile("data/RelLum.root");
  TTree *t_rlum = (TTree*)f_rlum->Get("T");
  int runnumber, spin_pattern;
  t_rlum->SetBranchAddress("Runnumber", &runnumber);
  t_rlum->SetBranchAddress("SpinPattern", &spin_pattern);

  TFile *f_yield = new TFile("data/yield-pattern.root", "RECREATE");
  TTree *t_yield = new TTree("T", "Yield");
  int crossing, ipt;
  double yield;
  t_yield->Branch("Runnumber", &runnumber, "Runnumber/I");
  t_yield->Branch("SpinPattern", &spin_pattern, "SpinPattern/I");
  t_yield->Branch("crossing", &crossing, "crossing/I");
  t_yield->Branch("ipt", &ipt, "ipt/I");
  t_yield->Branch("Yield", &yield, "Yield/D");

  for(int ien=0; ien<t_rlum->GetEntries(); ien++)
  {
    t_rlum->GetEntry(ien);

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/16669/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() )
    {
      cout << "Cannot open file for runnumber = " << runnumber << endl;
      delete f;
      continue;
    }

    int ical = 0;  // Use Sasha's calibration
    int isolated = 1;  // isolated photons
    int checkmap = 1;  // Use DC deadmap
    for(int icr=0; icr<2; icr++)
      for(int ib=0; ib<60; ib++)
      {
        crossing = icr + 2*ib;  // crossing ID in spin database
        int ih = isolated + 6*icr + 6*2*ib + 6*2*60*checkmap + 6*2*60*2*ical;
        TH1 *h_pt = (TH1*)f->Get(Form("h_photon_bunch_%d",ih));
        for(ipt=0; ipt<npT_pol; ipt++)
        {
          yield = h_pt->GetBinContent(ipt+1);
          t_yield->Fill();
        } // ipt
      } // icr, ib

    delete f;
  } // ien

  f_yield->cd();
  t_yield->Write();
  f_yield->Close();
  delete f_rlum;

  return 0;
}
