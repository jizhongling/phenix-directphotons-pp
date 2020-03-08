// To compile: g++ -Wall -o testIsoPhotonALL testIsoPhotonALL.cc `root-config --cflags --libs`
#include <iostream>
#include <sstream>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

using namespace std;

bool valid_bunch(int bunch, int spin_pol[])
{
  for(int i=0; i<4; i++)
    if( abs(spin_pol[(bunch-i+120)%120]) != 1 )
      return false;

  return true;
}

int main()
{
  const int taxi = 15811;
  const int npT_pol = 15;

  TFile *f_rlum = new TFile("data/RelLum.root");
  TTree *t_rlum = (TTree*)f_rlum->Get("T");
  const int nruns = t_rlum->GetEntries();
  int runnumber, fillnumber, spin_pattern;
  unsigned long long runevents, fillevents;
  double pol[2], epol[2], rate_bunch[120], erate_bunch[120];
  int spin_pol[120];
  double rlum[2], erlum[2];
  t_rlum->SetBranchAddress("Runnumber", &runnumber);
  t_rlum->SetBranchAddress("Fillnumber", &fillnumber);
  t_rlum->SetBranchAddress("SpinPattern", &spin_pattern);
  t_rlum->SetBranchAddress("RunEvents", &runevents);
  t_rlum->SetBranchAddress("FillEvents", &fillevents);
  t_rlum->SetBranchAddress("Pol", pol);
  t_rlum->SetBranchAddress("ePol", epol);
  t_rlum->SetBranchAddress("Pattern", spin_pol);
  t_rlum->SetBranchAddress("Rate", rate_bunch);
  t_rlum->SetBranchAddress("eRate", erate_bunch);

  TFile *f_ken2 = new TFile("data/YieldKEN2-isophoton.root");
  TTree *t_ken2 = (TTree*)f_ken2->Get("t1");
  const int nken2 = t_ken2->GetEntries();
  int ken2_ipt, ken2_part;
  double ken2_value;
  double *ken2 = new double[nken2];
  t_ken2->SetBranchAddress("ipt", &ken2_ipt);
  t_ken2->SetBranchAddress("part", &ken2_part);
  t_ken2->SetBranchAddress("value", &ken2_value);
  for(int ien=0; ien<nken2; ien++)
  {
    t_ken2->GetEntry(ien);
    int index = ken2_part + 6*3*2*2*2*ken2_ipt; 
    ken2[index] = ken2_value;
  }

  TFile *f_all = new TFile("data/isophoton-test.root", "RECREATE");
  TTree *t_all = new TTree("t1", "t1 tree");
  int ALLtype;
  double runno, ALL, eALL;
  t_all->Branch("ipt", &runnumber, "ipt/I");
  t_all->Branch("part", &ALLtype, "part/I");
  t_all->Branch("xpt", &runno, "xpt/D");
  t_all->Branch("value", &ALL, "value/D");
  t_all->Branch("error", &eALL, "error/D");
  t_all->Branch("errorlow", &eALL, "errorlow/D");
  t_all->Branch("errorhigh", &eALL, "errorhigh/D");

  for(int ien=0; ien<nruns; ien++)
  {
    t_rlum->GetEntry(ien);
    if( runevents < fillevents/3*2 )
      continue;

    runno = (double)runnumber;

    /* Calculate relative luminosities under this random spin patterns */
    double rate[2][2] = {}, e2rate[2][2] = {};  // icr, ipol
    for(int ib=0; ib<120; ib++)
      //if(abs(spin_pol[ib]) == 1)
      if(valid_bunch(ib,spin_pol))
      {
        int ipol = spin_pol[ib] > 0 ? 1 : 0;
        rate[ib%2][ipol] += rate_bunch[ib];
        e2rate[ib%2][ipol] += erate_bunch[ib]*erate_bunch[ib];
      }
    for(int icr=0; icr<2; icr++)
    {
      rlum[icr] = rate[icr][1]/rate[icr][0];
      erlum[icr] = rlum[icr]*sqrt(e2rate[icr][1]/rate[icr][1]/rate[icr][1] +
          e2rate[icr][0]/rate[icr][0]/rate[icr][0]);
    }

    /* Calculate ALL for each run */
    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/%d/data/PhotonHistos-%d.root",taxi,runnumber));
    if( f->IsZombie() )
    {
      cout << "Cannot open file for runnumber = " << runnumber << endl;
      delete f;
      continue;
    }

    int ical = 0;
    int checkmap = 1;
    int beam = 2;

    double nphoton[6][2][2][npT_pol] = {};  // imul, icr, ipol, ipt
    for(int imul=0; imul<6; imul++)
      for(int icr=0; icr<2; icr++)
        for(int ib=0; ib<60; ib++)
          //if(abs(spin_pol[icr + 2*ib]) == 1)
          if(valid_bunch(icr + 2*ib,spin_pol))
          {
            int ipol = spin_pol[icr + 2*ib] > 0 ? 1 : 0;
            int ih = imul + 6*icr + 6*2*ib + 6*2*60*checkmap + 6*2*60*2*ical;
            TH1 *h_photon_bunch = (TH1*)f->Get(Form("h_photon_bunch_%d",ih));
            for(int ipt=0; ipt<npT_pol; ipt++)
              nphoton[imul][icr][ipol][ipt] += h_photon_bunch->GetBinContent(ipt+1);
          } // icr, ib, imul
    for(int icr=0; icr<2; icr++)
      for(int ipol=0; ipol<2; ipol++)
        for(int ipt=0; ipt<npT_pol; ipt++)
          nphoton[0][icr][ipol][ipt] += nphoton[1][icr][ipol][ipt];

    for(int icr=0; icr<2; icr++)
    {
      double pbeam = pol[0]*pol[1];
      double e2pbeam = epol[0]*epol[0]/pol[0]/pol[0] + epol[1]*epol[1]/pol[1]/pol[1];
      double r = rlum[icr];
      double er = erlum[icr];

      for(int imul=0; imul<6; imul++)
        for(int ipt=0; ipt<npT_pol; ipt++)
        {
          double npp = nphoton[imul][icr][1][ipt];
          double npm = nphoton[imul][icr][0][ipt];

          int index = imul + 6*beam + 6*3*icr + 6*3*2*checkmap + 6*3*2*2*ical + 6*3*2*2*2*ipt;
          double k2 = ken2[index];
          ALL = 1./pbeam*(npp - r*npm)/(npp + r*npm);
          eALL = sqrt(pow(2*r*npp*npm/pbeam,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
              + e2pbeam*ALL*ALL);

          if( TMath::Finite(ALL+eALL) && eALL > 0. )
          {
            ALLtype = imul + 6*beam + 6*3*icr + 6*3*2*spin_pattern + 6*3*2*4*ipt;
            t_all->Fill();
          }
        } // imul, ipt
    } // icr

    delete f;
  } // ien

  f_all->cd();
  t_all->Write();
  f_all->Close();
  f_rlum->Close();
  f_ken2->Close();
  delete[] ken2;

  return 0;
}
