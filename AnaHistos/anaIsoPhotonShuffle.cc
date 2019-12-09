// To compile: g++ -Wall -o anaIsoPhotonShuffle anaIsoPhotonShuffle.cc `root-config --cflags --libs`
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TMath.h>
#include <TF1.h>
#include <TRandom3.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>

using namespace std;

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    cout << "Usage: " << argv[0] << " <Process>" << endl;
    exit(1);
  }

  int process;
  stringstream ss;
  ss << argv[1];
  ss >> process;

  /* pT bins for ALL */
  const int npT_pol = 15;
  const double pTbin_pol[npT_pol+1] = { 2.0,
    2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0,
    10.0, 12.0, 15.0, 20.0, 30.0 };

  TFile *f_rlum = new TFile("data/RelLum.root");
  TTree *t_rlum = (TTree*)f_rlum->Get("T");
  int runnumber, fillnumber, lastfill = 0;
  double pol[2], epol[2], rate_bunch[120], erate_bunch[120];
  double pattern[120], rlum[2], erlum[2];
  t_rlum->SetBranchAddress("Runnumber", &runnumber);
  t_rlum->SetBranchAddress("Fillnumber", &fillnumber);
  t_rlum->SetBranchAddress("Pol", pol); t_rlum->SetBranchAddress("ePol", epol);
  t_rlum->SetBranchAddress("Rate", rate_bunch);
  t_rlum->SetBranchAddress("eRate", erate_bunch);

  TFile *f_all = new TFile(Form("histos/isophoton-shuffle-%d.root",process), "RECREATE");
  TTree *t_all = new TTree("T", "ALL bunch shuffling");
  int type;
  double ALLmean, eALLmean;
  t_all->Branch("Process", &process, "Process/I");
  t_all->Branch("Type", &type, "Type/I");
  t_all->Branch("ALL", &ALLmean, "ALL/D");
  t_all->Branch("eALL", &eALLmean, "eALL/D");

  const int ngr_asym = 5*2*npT_pol;
  TGraphErrors *gr_asym[ngr_asym];
  int igr_asym[ngr_asym] = {};
  for(int ig=0; ig<ngr_asym; ig++)
    gr_asym[ig] = new TGraphErrors(800);

  TRandom3 *rnd = new TRandom3();

  TF1 *fn_mean = new TF1("fn_mean", "pol0");

  //for(int ien=0; ien<t_rlum->GetEntries(); ien++)
  for(int ien=0; ien<10; ien++)
  {
    t_rlum->GetEntry(ien);

    /* Generate spin patterns for each fill */
    if( fillnumber != lastfill )
    {
      double rate[2][2] = {}, e2rate[2][2] = {};
      rnd->SetSeed(ien + 1000*process + 1);
      for(int ib=0; ib<120; ib++)
      {
        int ipol = rnd->Integer(2);
        pattern[ib] = ipol;
        rate[ib%2][ipol] += rate_bunch[ib];
        e2rate[ib%2][ipol] += erate_bunch[ib]*erate_bunch[ib];
      }

      for(int icr=0; icr<2; icr++)
      {
        rlum[icr] = rate[icr][0]/rate[icr][1];
        erlum[icr] = rlum[icr]*sqrt(e2rate[icr][0]/rate[icr][0]/rate[icr][0] +
            e2rate[icr][1]/rate[icr][1]/rate[icr][1]);
      }

      lastfill = fillnumber;
    }

    /* Calculate ALL for each run */
    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15717/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() )
    {
      cout << "Cannot open file for runnumber = " << runnumber << endl;
      continue;
    }

    double nphoton[2][2][npT_pol] = {};  // icr, ipol, ipt
    double npion[2][2][2][2][npT_pol] = {};  // pttype(pt|2pt), ibg, icr, ipol, ipt

    int checkmap = 1;

    for(int icr=0; icr<2; icr++)
      for(int ib=0; ib<60; ib++)
      {
        int ipol = pattern[ib*2];
        int ih = icr + 2*ib + 2*60*checkmap;
        TH1 *h_photon = (TH1*)f->Get(Form("h_1photon_bunch_%d",ih));
        for(int ipt=0; ipt<npT_pol; ipt++)
        {
          int ptbin_first = h_photon->GetXaxis()->FindBin(pTbin_pol[ipt]);
          int ptbin_last = h_photon->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;
          nphoton[icr][ipol][ipt] += h_photon->Integral(ptbin_first,ptbin_last);
        } // ipt

        for(int pttype=0; pttype<2; pttype++)
        {
          const char *ptname = pttype ? "2pt" : "";
          int ih = icr + 2*ib + 2*60*checkmap;
          TH2 *h2_pion = (TH2*)f->Get(Form("h2_2photon%s_bunch_%d",ptname,ih));
          for(int ipt=0; ipt<npT_pol; ipt++)
          {
            int ptbin_first = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt]);
            int ptbin_last = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;
            npion[pttype][0][icr][ipol][ipt] += h2_pion->Integral(ptbin_first,ptbin_last, 113,162);
            npion[pttype][1][icr][ipol][ipt] += h2_pion->Integral(ptbin_first,ptbin_last, 48,97) +
              h2_pion->Integral(ptbin_first,ptbin_last, 178,227);
          } // ipt
        } // pttype
      } // icr, ib

    for(int icr=0; icr<2; icr++)
    {
      double pbeam = pol[0]*pol[1];
      double e2pbeam = epol[0]*epol[0]/pol[0]/pol[0] + epol[1]*epol[1]/pol[1]/pol[1];
      double r = rlum[icr];
      double er = erlum[icr];

      for(int ipt=0; ipt<npT_pol; ipt++)
      {
        for(int pttype=0; pttype<2; pttype++)
          for(int ibg=0; ibg<2; ibg++)
          {
            double npp = npion[pttype][ibg][icr][1][ipt];
            double npm = npion[pttype][ibg][icr][0][ipt];

            double k2 = 1.;
            double ALL = 1./pbeam*(npp - r*npm)/(npp + r*npm);
            double eALL = sqrt(pow(2*r*npp*npm/pbeam,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
                + e2pbeam*ALL*ALL);

            if( TMath::Finite(ALL+eALL) && eALL > 0. )
            {
              int imul = 1 + ibg + 2*pttype;
              int ig = imul + 5*icr + 5*2*ipt;
              gr_asym[ig]->SetPoint(igr_asym[ig], (double)runnumber, ALL);
              gr_asym[ig]->SetPointError(igr_asym[ig], 0., eALL);
              igr_asym[ig]++;
            }
          } // pttype, ibg

        double npp = nphoton[icr][1][ipt];
        double npm = nphoton[icr][0][ipt];

        double k2 = 1.;
        double ALL = 1./pbeam*(npp - r*npm)/(npp + r*npm);
        double eALL = sqrt(pow(2*r*npp*npm/pbeam,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
            + e2pbeam*ALL*ALL);

        if( TMath::Finite(ALL+eALL) && eALL > 0. )
        {
          int imul = 0;
          int ig = imul + 5*icr + 5*2*ipt;
          gr_asym[ig]->SetPoint(igr_asym[ig], (double)runnumber, ALL);
          gr_asym[ig]->SetPointError(igr_asym[ig], 0., eALL);
          igr_asym[ig]++;
        }
      } // ipt
    } // icr

    delete f;
  } // ien

  /* Fit ALL for different regions and crossings */
  for(int ig=0; ig<ngr_asym; ig++)
  {
    type = ig;
    gr_asym[ig]->Set(igr_asym[ig]);

    gr_asym[ig]->Fit(fn_mean, "Q");
    ALLmean = fn_mean->GetParameter(0);
    eALLmean = fn_mean->GetParError(0);
    if( TMath::Finite(ALLmean+eALLmean) )
      t_all->Fill();
  }

  f_all->cd();
  t_all->Write();
  f_all->Close();

  return 0;
}
