// To compile: g++ -Wall -o anaIsoPhotonShuffle anaIsoPhotonShuffle.cc `root-config --cflags --libs`
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TMath.h>
#include <TFitResult.h>
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

  const int npT_pol = 15;

  TFile *f_rlum = new TFile("data/RelLum.root");
  TTree *t_rlum = (TTree*)f_rlum->Get("T");
  const int nruns = t_rlum->GetEntries();
  int runnumber, fillnumber, lastfill = 0;
  double pol[2], epol[2], rate_bunch[120], erate_bunch[120];
  int spin_pol[120], spin_rnd[120];
  double rlum[2], erlum[2];
  t_rlum->SetBranchAddress("Runnumber", &runnumber);
  t_rlum->SetBranchAddress("Fillnumber", &fillnumber);
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
    int index = ken2_part + 5*3*2*2*ken2_ipt; 
    ken2[index] = ken2_value;
  }

  TFile *f_all = new TFile(Form("histos/isophoton-shuffle-%d.root",process), "RECREATE");
  TTree *t_all = new TTree("T", "ALL bunch shuffling");
  int ALLtype, ALLndf;
  double ALLchi2, ALLmean, eALLmean;
  t_all->Branch("Process", &process, "Process/I");
  t_all->Branch("Type", &ALLtype, "Type/I");
  t_all->Branch("NDF", &ALLndf, "NDF/I");
  t_all->Branch("Chi2", &ALLchi2, "Chi2/D");
  t_all->Branch("ALL", &ALLmean, "ALL/D");
  t_all->Branch("eALL", &eALLmean, "eALL/D");

  const int ngr_asym = 5*2*npT_pol;
  TGraphErrors *gr_asym[ngr_asym];
  int igr_asym[ngr_asym] = {};
  for(int ig=0; ig<ngr_asym; ig++)
    gr_asym[ig] = new TGraphErrors(nruns);

  TRandom3 *rnd = new TRandom3();

  for(int ien=0; ien<nruns; ien++)
  {
    t_rlum->GetEntry(ien);

    /* Generate random spin patterns for each fill */
    if(fillnumber != lastfill)
    {
      rnd->SetSeed(ien + 1000*process + 1); 
      for(int ib=0; ib<120; ib++)
        if(abs(spin_pol[ib]) == 1)
          spin_rnd[ib] = rnd->Integer(2);
      lastfill = fillnumber;
    }

    /* Calculate relative luminosities under this random spin patterns */
    double rate[2][2] = {}, e2rate[2][2] = {};  // icr, ipol
    for(int ib=0; ib<120; ib++)
      if(abs(spin_pol[ib]) == 1)
      {
        int ipol = spin_rnd[ib];
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
    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15763/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() )
    {
      cout << "Cannot open file for runnumber = " << runnumber << endl;
      delete f;
      continue;
    }

    double nphoton[2][2][npT_pol] = {};  // icr, ipol, ipt
    double npion[2][2][2][2][npT_pol] = {};  // pttype(pt|2pt), ibg, icr, ipol, ipt

    int checkmap = 1;
    int beam = 2;

    for(int icr=0; icr<2; icr++)
      for(int ib=0; ib<60; ib++)
        if(abs(spin_pol[icr + 2*ib]) == 1)
        {
          int ipol = spin_rnd[icr + 2*ib];
          int ih = icr + 2*ib + 2*60*checkmap;
          // TODO: future TAXI will have different names
          TH1 *h_photon = (TH1*)f->Get(Form("h_1photon_bunch_%d",ih));
          for(int ipt=0; ipt<npT_pol; ipt++)
            nphoton[icr][ipol][ipt] += h_photon->GetBinContent(ipt+1);

          for(int pttype=0; pttype<2; pttype++)
          {
            const char *ptname = pttype ? "2pt" : "";
            int ih = icr + 2*ib + 2*60*checkmap;
            // TODO: future TAXI will have different names
            TH2 *h2_pion = (TH2*)f->Get(Form("h2_2photon%s_bunch%d",ptname,ih));
            for(int ipt=0; ipt<npT_pol; ipt++)
            {
              npion[pttype][0][icr][ipol][ipt] += h2_pion->GetBinContent(ipt+1,3);
              npion[pttype][1][icr][ipol][ipt] += h2_pion->GetBinContent(ipt+1,1) +
                h2_pion->GetBinContent(ipt+1,5);
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

            int imul = 1 + ibg + 2*pttype;
            int index = imul + 5*beam + 5*3*icr + 5*3*2*checkmap + 5*3*2*2*ipt;
            double k2 = ken2[index];
            double ALL = 1./pbeam*(npp - r*npm)/(npp + r*npm);
            double eALL = sqrt(pow(2*r*npp*npm/pbeam,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
                + e2pbeam*ALL*ALL);

            if( TMath::Finite(ALL+eALL) && eALL > 0. )
            {
              int ig = imul + 5*icr + 5*2*ipt;
              gr_asym[ig]->SetPoint(igr_asym[ig], (double)runnumber, ALL);
              gr_asym[ig]->SetPointError(igr_asym[ig], 0., eALL);
              igr_asym[ig]++;
            }
          } // pttype, ibg

        double npp = nphoton[icr][1][ipt];
        double npm = nphoton[icr][0][ipt];

        int imul = 0;
        int index = imul + 5*beam + 5*3*icr + 5*3*2*checkmap + 5*3*2*2*ipt;
        double k2 = ken2[index];
        double ALL = 1./pbeam*(npp - r*npm)/(npp + r*npm);
        double eALL = sqrt(pow(2*r*npp*npm/pbeam,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
            + e2pbeam*ALL*ALL);

        if( TMath::Finite(ALL+eALL) && eALL > 0. )
        {
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
    ALLtype = ig;
    gr_asym[ig]->Set(igr_asym[ig]);

    TFitResultPtr r_asym = gr_asym[ig]->Fit("pol0", "QS");
    ALLndf = r_asym->Ndf();
    if(ALLndf > 0)
    {
      ALLchi2 = r_asym->Chi2()/ALLndf;
      ALLmean = r_asym->Value(0);
      eALLmean = r_asym->ParError(0);
      if( TMath::Finite(ALLmean+eALLmean) )
        t_all->Fill();
    }

    delete gr_asym[ig];
  }

  f_all->cd();
  t_all->Write();
  f_all->Close();
  f_rlum->Close();
  f_ken2->Close();
  delete[] ken2;
  delete rnd;

  return 0;
}
