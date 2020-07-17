// To compile: g++ -Wall -o calcRelLum calcRelLum.cc -I$OFFLINE_MAIN/include -L$OFFLINE_MAIN/lib -m32 -luspin -lodbc -lgsl -lgslcblas -lm `root-config --cflags --libs`
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/foreach.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>

#include <SpinDBOutput.hh>
#include <SpinDBContent.hh>

#include "CommonFunc.h"

using namespace std;

struct Pol
{
  double pol[2];
  double epol[2];
};

double rate_f(double x, void *params)
{
  double *p = (double*)params;
  double kn = p[0];
  double ks = p[1];
  double rate_obs = p[2];

  double rate_true = 1 - exp(-x*(1+kn)) - exp(-x*(1+ks)) + exp(-x*(1+kn+ks));

  return rate_true - rate_obs;
}

double rate_df(double x, void *params)
{
  double *p = (double*)params;
  double kn = p[0];
  double ks = p[1];

  double df = (1+kn)*exp(-x*(1+kn)) + (1+ks)*exp(-x*(1+ks)) - (1+kn+ks)*exp(-x*(1+kn+ks));

  return df;
}

void rate_fdf(double x, void *params, double *y, double *dy)
{
  double *p = (double*)params;
  double kn = p[0];
  double ks = p[1];
  double rate_obs = p[2];

  double expkn = exp(-x*(1+kn));
  double expks = exp(-x*(1+ks));
  double expknks = exp(-x*(1+kn+ks));

  double rate_true = 1 - expkn - expks + expknks;
  double df = (1+kn)*expkn + (1+ks)*expks - (1+kn+ks)*expknks;

  *y = rate_true - rate_obs;
  *dy = df;
}

// Input: KFactor[2] = {kn, ks}
// Input: rate_obs[2] = {rate_obs_vtx, rate_obs_novtx}
// Input: erate_obs[2] = {errors for rate_obs}
// Output: rate_true[3] = {rate_true_vtx, rate_true_novtx, rate_true_res}
// Output: erate_true[3] = {errors for rate_ture}
void rate_corr(double KFactor[], double rate_obs[], double erate_obs[],
    double rate_true[], double erate_true[])
{
  const int max_iter = 100;
  const double epsrel = 1e-6;

  const gsl_root_fdfsolver_type *solver_type = gsl_root_fdfsolver_newton;
  gsl_root_fdfsolver *solver = gsl_root_fdfsolver_alloc(solver_type);

  double params[3] = {KFactor[0], KFactor[1]};
  gsl_function_fdf FDF;
  FDF.f = &rate_f;
  FDF.df = &rate_df;
  FDF.fdf = &rate_fdf;
  FDF.params = params;

  double erel_obs[2], erel_true[3];

  for(int i=0; i<3; i++)
  {
    if(i < 2)
      params[2] = rate_obs[i];
    else
      params[2] = rate_obs[0]*(rate_obs[0]/rate_obs[1])/(rate_true[0]/rate_true[1]);

    int status;
    int iter = 0;
    double x0, x = params[2];
    gsl_root_fdfsolver_set(solver, &FDF, x);

    do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate(solver);
      x0 = x;
      x = gsl_root_fdfsolver_root(solver);
      status = gsl_root_test_delta(x, x0, 0., epsrel);
    }
    while(status == GSL_CONTINUE && iter < max_iter);

    if(status != GSL_SUCCESS)
    {
      cerr << "Not converged!" << endl;
      gsl_root_fdfsolver_free(solver);
      exit(1);
    }

    rate_true[i] = x;
    if(i < 2)
    {
      erel_obs[i] = erate_obs[i]/params[2];
      erel_true[i] = erel_obs[i];
    }
    else
    {
      erel_true[i] = sqrt(4*erel_obs[0]*erel_obs[0] + erel_obs[1]*erel_obs[1] +
          erel_true[0]*erel_true[0] + erel_true[1]*erel_true[1]);
    }
    erel_true[i] *= params[2]/x/rate_df(x,params);
    erate_true[i] = erel_true[i]*x;
  }

  gsl_root_fdfsolver_free(solver);
  return;
}

void get_gl1p(SpinDBOutput &spin_out, SpinDBContent &spin_cont,
    int runnumber, int &fillnumber, int &lastfill,
    double pol[], double epol[], double rlum[][2], double erlum[][2],
    int spin_pol[], double count_bunch[], double ecount_bunch[],
    double count_fill[][2][2], double e2count_fill[][2][2])
{
  /* Retrieve entry from Spin DB and get fill number */
  int qa_level = spin_out.GetDefaultQA(runnumber);
  spin_out.StoreDBContent(runnumber, runnumber, qa_level);
  spin_out.GetDBContentStore(spin_cont, runnumber);
  fillnumber = spin_cont.GetFillNumber();
  if(fillnumber != lastfill)
  {
    lastfill = fillnumber;
    for(int beam=0; beam<3; beam++)
      for(int evenodd=0; evenodd<2; evenodd++)
        for(int ipol=0; ipol<2; ipol++)
        {
          count_fill[beam][evenodd][ipol] = 0.;
          e2count_fill[beam][evenodd][ipol] = 0.;
        }
  }

  if( spin_out.CheckRunRow(runnumber,qa_level) != 1 )
  {
    cerr << "Wrong spinDB info!" << endl;
    exit(1);
  }

  double pbsyst;
  double pysyst;
  spin_cont.GetPolarizationBlue(1, pol[0], epol[0], pbsyst);
  spin_cont.GetPolarizationYellow(1, pol[1], epol[1], pysyst);

  long long bbc30_gl1p[3][2][2] = {};

  for(int ib=0; ib<120; ib++)
  {
    count_bunch[ib] = 0.;
    ecount_bunch[ib] = 0.;
    int pb = spin_cont.GetSpinPatternBlue(ib);
    int py = spin_cont.GetSpinPatternYellow(ib);
    int pattern[3] = {pb, py, pb*py};
    for(int beam=0; beam<3; beam++)
      if( abs(pattern[beam]) == 1 )
      {
        int ipol = pattern[beam] > 0 ? 1 : 0;
        long long count = spin_cont.GetScaler(0, ib);
        bbc30_gl1p[beam][ib%2][ipol] += count;
        if(beam == 2)
        {
          if( abs(pattern[2]) == 1 )
            spin_pol[ib] = pattern[2];
          else
            spin_pol[ib] = 0;
          count_bunch[ib] = (double)count;
          ecount_bunch[ib] = sqrt(count_bunch[ib]);
        }
      }
  }

  for(int beam=0; beam<3; beam++)
    for(int evenodd=0; evenodd<2; evenodd++)
    {
      rlum[beam][evenodd] = (double)bbc30_gl1p[beam][evenodd][1]/bbc30_gl1p[beam][evenodd][0];
      erlum[beam][evenodd] = rlum[beam][evenodd]*sqrt(1./bbc30_gl1p[beam][evenodd][1] + 1./bbc30_gl1p[beam][evenodd][0]);
      for(int ipol=0; ipol<2; ipol++)
      {
        count_fill[beam][evenodd][ipol] += bbc30_gl1p[beam][evenodd][ipol];
        e2count_fill[beam][evenodd][ipol] += bbc30_gl1p[beam][evenodd][ipol];
      }
    }

  return;
}

void ReadKFactor(double KFactor[][2])
{
  TFile *f_k = new TFile("/phenix/plhf/zji/sources/offline/AnalysisTrain/Run13_Pi0Ana_YIS/ResultKFactor.root");

  for(int i=0; i<5; i++)
  {
    string pattern;
    if(i==0) pattern = "SOOSSOO";
    else if(i==1) pattern = "OSSOOSS";
    else if(i==2) pattern = "SSOO";
    else if(i==3) pattern = "OOSS";
    else if(i==4) pattern = "ALL";

    TGraphErrors *gr_bbc_north = (TGraphErrors*)f_k->Get(Form("Rate_Ratio_Corr_BBC_North_%s_0", pattern.c_str()));
    TGraphErrors *gr_bbc_south = (TGraphErrors*)f_k->Get(Form("Rate_Ratio_Corr_BBC_South_%s_0", pattern.c_str()));

    TF1 *fr_bbc_north = gr_bbc_north->GetFunction(Form("BBC_NORTH_CORR_%s", pattern.c_str()));
    TF1 *fr_bbc_south = gr_bbc_south->GetFunction(Form("BBC_SOUTH_CORR_%s", pattern.c_str()));

    KFactor[i][0] = fr_bbc_north->GetParameter(0);
    KFactor[i][1] = fr_bbc_south->GetParameter(0);
  }

  delete f_k;
  return;
}

int main()
{
  int runnumber, fillnumber, lastfill = 0;
  int spin_pattern = 0;
  unsigned long long runevents, fillevents;

  Pol pol_fill;
  map<int,Pol> pol_info;
  ifstream fin("data/pol-fill.txt");
  while(fin >> fillnumber >> pol_fill.pol[0] >> pol_fill.epol[0] >> pol_fill.pol[1] >> pol_fill.epol[1])
  {
    for(int i=0; i<2; i++)
    {
      pol_fill.pol[i] /= 100.;
      pol_fill.epol[i] /= 100.;
    }
    pol_info.insert( make_pair(fillnumber, pol_fill) );
  }
  fin.close();

  TFile *f_rlum = new TFile("data/RelLum.root", "RECREATE");

  TTree *t_rlum = new TTree("T", "Relative luminosity");
  int spin_pol[120];
  double pol[2], epol[2], rlum[3][2], erlum[3][2],
         rlum_fill[3][2], erlum_fill[3][2],
         count_bunch[120], ecount_bunch[120];
  t_rlum->Branch("Runnumber", &runnumber, "Runnumber/I");
  t_rlum->Branch("Fillnumber", &fillnumber, "Fillnumber/I");
  t_rlum->Branch("SpinPattern", &spin_pattern, "SpinPattern/I");
  t_rlum->Branch("RunEvents", &runevents, "RunEvents/l");
  t_rlum->Branch("FillEvents", &fillevents, "FillEvents/l");
  t_rlum->Branch("Pol", pol, "Pol[2]/D");
  t_rlum->Branch("ePol", epol, "ePol[2]/D");
  t_rlum->Branch("PolFill", pol_fill.pol, "PolFill[2]/D");
  t_rlum->Branch("ePolFill", pol_fill.epol, "ePolFill[2]/D");
  t_rlum->Branch("BlueRelLum", rlum[0], "BlueRelLum[2]/D");
  t_rlum->Branch("eBlueRelLum", erlum[0], "eBlueRelLum[2]/D");
  t_rlum->Branch("YellowRelLum", rlum[1], "YellowRelLum[2]/D");
  t_rlum->Branch("eYellowRelLum", erlum[1], "eYellowRelLum[2]/D");
  t_rlum->Branch("RelLum", rlum[2], "RelLum[2]/D");
  t_rlum->Branch("eRelLum", erlum[2], "eRelLum[2]/D");
  t_rlum->Branch("BlueRelLumFill", rlum_fill[0], "BlueRelLumFill[2]/D");
  t_rlum->Branch("eBlueRelLumFill", erlum_fill[0], "eBlueRelLumFill[2]/D");
  t_rlum->Branch("YellowRelLumFill", rlum_fill[1], "YellowRelLumFill[2]/D");
  t_rlum->Branch("eYellowRelLumFill", erlum_fill[1], "eYellowRelLumFill[2]/D");
  t_rlum->Branch("RelLumFill", rlum_fill[2], "RelLumFill[2]/D");
  t_rlum->Branch("eRelLumFill", erlum_fill[2], "eRelLumFill[2]/D");
  t_rlum->Branch("Pattern", spin_pol, "Pattern[120]/I");
  t_rlum->Branch("Count", count_bunch, "Count[120]/D");
  t_rlum->Branch("eCount", ecount_bunch, "eCount[120]/D");

  TTree *t_rlum_check = new TTree("T1", "Relative luminosity check");
  bool runqa;
  double gl1p[3][2], egl1p[3][2], ss_raw[2], ess_raw[2],
         ss_pile[2], ess_pile[2], ss_res[2], ess_res[2];
  t_rlum_check->Branch("Runnumber", &runnumber, "Runnumber/I");
  t_rlum_check->Branch("RunQA", &runqa, "RunQA/O");
  t_rlum_check->Branch("GL1p", gl1p[2], "Gl1p[2]/D");
  t_rlum_check->Branch("eGL1p", egl1p[2], "eGl1p[2]/D");
  t_rlum_check->Branch("SS_Uncorr", ss_raw, "SS_Uncorr[2]/D");
  t_rlum_check->Branch("eSS_Uncorr", ess_raw, "eSS_Uncorr[2]/D");
  t_rlum_check->Branch("SS_Pileup", ss_pile, "SS_Pileup[2]/D");
  t_rlum_check->Branch("eSS_Pileup", ess_pile, "eSS_Pileup[2]/D");
  t_rlum_check->Branch("SS_Residual", ss_res, "SS_Residual[2]/D");
  t_rlum_check->Branch("eSS_Residual", ess_res, "eSS_Residual[2]/D");

  typedef map<int,int> map_int_t;
  map_int_t runnoInseok;  // <runnumber, spin_pattern>
  ifstream finInseok("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist-ALL.txt");
  while( finInseok >> runnumber )
  {
    if(runnumber == 0) { spin_pattern++; continue; }
    if(spin_pattern > 3) break;
    runnoInseok.insert( make_pair(runnumber, spin_pattern) );
  }
  finInseok.close();

  SpinDBOutput spin_out;
  SpinDBContent spin_cont;

  /* Initialize object to access spin DB */
  spin_out.Initialize();
  spin_out.SetUserName("phnxrc");
  spin_out.SetTableName("spin");

  TFile *f_bbc = new TFile("data/BBC_Run13_PING.root");
  TTree *t_bbc = (TTree*)f_bbc->Get("T");
  int pb[120], py[120], clock[120], bbc30[120], bbcnovtx[120];
  t_bbc->SetBranchAddress("Runnumber", &runnumber);
  t_bbc->SetBranchAddress("SpinPatBlue", pb);
  t_bbc->SetBranchAddress("SpinPatYellow", py);
  t_bbc->SetBranchAddress("CLOCK", clock);
  t_bbc->SetBranchAddress("BBC_30cm", bbc30);
  t_bbc->SetBranchAddress("BBCnovtx", bbcnovtx);

  typedef map<int,unsigned long long> map_ulong_t;
  map_ulong_t daq_events, daq_runevents, daq_fillevents;  // <runnumber, nevents>
  TFile *f_daq = new TFile("data/clock-counts.root");
  TTree *t_daq = (TTree*)f_daq->Get("t1");
  unsigned long long runnoDaq, bbcnarrow;
  t_daq->SetBranchAddress("runnumber", &runnoDaq);
  t_daq->SetBranchAddress("bbcnarrow_live", &bbcnarrow);
  vector<int> runlistDaq;
  for(int ien=0; ien<t_daq->GetEntries(); ien++)
  {
    t_daq->GetEntry(ien);
    runnumber = (int)runnoDaq;
    daq_events.insert( make_pair(runnumber, bbcnarrow) );
    runlistDaq.push_back(runnumber);
  }
  sort(runlistDaq.begin(), runlistDaq.end());
  lastfill = 0;
  vector<int> unfilledrun;
  BOOST_FOREACH(const int &runnumber, runlistDaq)
  {
    /* Retrieve entry from Spin DB and get fill number */
    int qa_level = spin_out.GetDefaultQA(runnumber);
    spin_out.StoreDBContent(runnumber, runnumber, qa_level);
    spin_out.GetDBContentStore(spin_cont, runnumber);
    fillnumber = spin_cont.GetFillNumber();

    map_ulong_t::iterator it_run = daq_events.find(runnumber);
    runevents = it_run != daq_events.end() ? it_run->second : 0;
    if(fillnumber != lastfill)
    {
      BOOST_FOREACH(const int &unfill, unfilledrun)
        daq_fillevents.insert( make_pair(unfill, fillevents) );
      unfilledrun.clear();
      fillevents = runevents;
    }
    else
      fillevents += runevents;

    daq_runevents.insert( make_pair(runnumber, fillevents) );
    unfilledrun.push_back(runnumber);
    lastfill = fillnumber;
  }
  BOOST_FOREACH(const int &unfill, unfilledrun)
    daq_fillevents.insert( make_pair(unfill, fillevents) );

  double KFactor[2] = {0.2278, 0.2278};
  double KFactorInseok[5][2];
  ReadKFactor(KFactorInseok);

  vector<int> runnoDone;

  double count_fill[3][2][2] = {};
  double e2count_fill[3][2][2] = {};
  double gl1p_fill[3][2][2] = {};
  double e2gl1p_fill[3][2][2] = {};

  lastfill = 0;
  for(int ien=0; ien<t_bbc->GetEntries(); ien++)
  {
    t_bbc->GetEntry(ien);

    map_int_t::iterator it_Inseok = runnoInseok.find(runnumber);
    spin_pattern = it_Inseok != runnoInseok.end() ? it_Inseok->second : 4;

    get_gl1p(spin_out, spin_cont, runnumber, fillnumber, lastfill,
        pol, epol, gl1p, egl1p, spin_pol, count_bunch, ecount_bunch,
        gl1p_fill, e2gl1p_fill);

    map_ulong_t::iterator it_run = daq_runevents.find(runnumber);
    map_ulong_t::iterator it_fill = daq_fillevents.find(runnumber);
    runevents = it_run != daq_runevents.end() ? it_run->second : 0;
    fillevents = it_fill != daq_fillevents.end() ? it_fill->second : 0;

    double count_raw[2][2] = {};
    double e2count_raw[2][2] = {};
    double count_pile[2][2] = {};
    double e2count_pile[2][2] = {};
    double count_res[3][2][2] = {};
    double e2count_res[3][2][2] = {};

    for(int ib=0; ib<120; ib++)
    {
      int pattern[3] = {pb[ib], py[ib], pb[ib]*py[ib]};

      if( clock[ib] > 0 && bbc30[ib] > 0 && bbcnovtx[ib] > 0 )
      {
        double rate_vtx = (double)bbc30[ib]/clock[ib];
        double erate_vtx = sqrt((double)bbc30[ib])/clock[ib];
        double rate_novtx = (double)bbcnovtx[ib]/clock[ib];
        double erate_novtx = sqrt((double)bbcnovtx[ib])/clock[ib];

        double rate_obs[2] = {rate_vtx, rate_novtx};
        double erate_obs[2] = {erate_vtx, erate_novtx};
        double rate_true[3], erate_true[3];
        rate_corr(KFactor, rate_obs, erate_obs, rate_true, erate_true);

        for(int beam=0; beam<3; beam++)
          if( abs(pattern[beam]) == 1 )
          {
            int ipol = pattern[beam] > 0 ? 1 : 0;
            if(beam == 2)
            {
              count_raw[ib%2][ipol] += rate_obs[0]*clock[ib];
              e2count_raw[ib%2][ipol] += pow2(erate_obs[0]*clock[ib]);
              count_pile[ib%2][ipol] += rate_true[0]*clock[ib];
              e2count_pile[ib%2][ipol] += pow2(erate_true[0]*clock[ib]);
              count_bunch[ib] = rate_true[2]*clock[ib];
              ecount_bunch[ib] = erate_true[2]*clock[ib];
            }
            count_res[beam][ib%2][ipol] += rate_true[2]*clock[ib];
            e2count_res[beam][ib%2][ipol] += pow2(erate_true[2]*clock[ib]);
          }
      }
      else
      {
        for(int beam=0; beam<3; beam++)
          if( abs(pattern[beam]) == 1 )
          {
            int ipol = pattern[beam] > 0 ? 1 : 0;
            count_res[beam][ib%2][ipol] += count_bunch[ib];
            e2count_res[beam][ib%2][ipol] += pow2(ecount_bunch[ib]);
          }
      }
    } // ib

    for(int evenodd=0; evenodd<2; evenodd++)
    {
      ss_raw[evenodd] = count_raw[evenodd][1]/count_raw[evenodd][0];
      ess_raw[evenodd] = ss_raw[evenodd]*sqrt(e2count_raw[evenodd][1]/pow2(count_raw[evenodd][1])
          + e2count_raw[evenodd][0]/pow2(count_raw[evenodd][0]));
      ss_pile[evenodd] = count_pile[evenodd][1]/count_pile[evenodd][0];
      ess_pile[evenodd] = ss_pile[evenodd]*sqrt(e2count_pile[evenodd][1]/pow2(count_pile[evenodd][1])
          + e2count_pile[evenodd][0]/pow2(count_pile[evenodd][0]));
      for(int beam=0; beam<3; beam++)
      {
        rlum[beam][evenodd] = count_res[beam][evenodd][1]/count_res[beam][evenodd][0];
        erlum[beam][evenodd] = rlum[beam][evenodd]*sqrt(e2count_res[beam][evenodd][1]/pow2(count_res[beam][evenodd][1])
            + e2count_res[beam][evenodd][0]/pow2(count_res[beam][evenodd][0]));
        for(int ipol=0; ipol<2; ipol++)
        {
          count_fill[beam][evenodd][ipol] += count_res[beam][evenodd][ipol];
          e2count_fill[beam][evenodd][ipol] += e2count_res[beam][evenodd][ipol];
        }
        rlum_fill[beam][evenodd] = count_fill[beam][evenodd][1]/count_fill[beam][evenodd][0];
        erlum_fill[beam][evenodd] = rlum_fill[beam][evenodd]*sqrt(e2count_fill[beam][evenodd][1]/pow2(count_fill[beam][evenodd][1])
            + e2count_fill[beam][evenodd][0]/pow2(count_fill[beam][evenodd][0]));
        if(beam == 2)
        {
          ss_res[evenodd] = rlum[beam][evenodd];
          ess_res[evenodd] = erlum[beam][evenodd];
        }
      }
    }

    runqa = false;
    if( it_Inseok != runnoInseok.end() )
    {
      map<int,Pol>::iterator it_pol = pol_info.find(fillnumber);
      if(it_pol == pol_info.end() )
        cout << "No pol info for fill " << fillnumber << endl;
      else
        for(int i=0; i<2; i++)
        {
          pol_fill.pol[i] = it_pol->second.pol[i];
          pol_fill.epol[i] = it_pol->second.epol[i];
        }

      if(runnumber < 386946)
        for(int beam=0; beam<3; beam++)
          for(int evenodd=0; evenodd<2; evenodd++)
          {
            rlum[beam][evenodd] = gl1p[beam][evenodd];
            erlum[beam][evenodd] = egl1p[beam][evenodd];
            rlum_fill[beam][evenodd] = gl1p_fill[beam][evenodd][1]/gl1p_fill[beam][evenodd][0];
            erlum_fill[beam][evenodd] = rlum_fill[beam][evenodd]*sqrt(e2gl1p_fill[beam][evenodd][1]/pow2(gl1p_fill[beam][evenodd][1])
                + e2gl1p_fill[beam][evenodd][0]/pow2(gl1p_fill[beam][evenodd][0]));
          }
      else
        runqa = true;

      t_rlum->Fill();
      runnoDone.push_back(runnumber);
    }

    t_rlum_check->Fill();
  } // ien

  lastfill = 0;
  BOOST_FOREACH(const map_int_t::value_type &runno, runnoInseok)
    if( find(runnoDone.begin(), runnoDone.end(), runno.first) == runnoDone.end() )
    {
      runnumber = runno.first;
      spin_pattern = runno.second;

      get_gl1p(spin_out, spin_cont, runnumber, fillnumber, lastfill,
          pol, epol, rlum, erlum, spin_pol, count_bunch, ecount_bunch,
          count_fill, e2count_fill);

      map_ulong_t::iterator it_run = daq_runevents.find(runnumber);
      map_ulong_t::iterator it_fill = daq_fillevents.find(runnumber);
      runevents = it_run != daq_runevents.end() ? it_run->second : 0;
      fillevents = it_fill != daq_fillevents.end() ? it_fill->second : 0;

      map<int,Pol>::iterator it_pol = pol_info.find(fillnumber);
      if(it_pol == pol_info.end() )
        cout << "No pol info for fill " << fillnumber << endl;
      else
        for(int i=0; i<2; i++)
        {
          pol_fill.pol[i] = it_pol->second.pol[i];
          pol_fill.epol[i] = it_pol->second.epol[i];
        }

      for(int beam=0; beam<3; beam++)
        for(int evenodd=0; evenodd<2; evenodd++)
        {
          rlum_fill[beam][evenodd] = count_fill[beam][evenodd][1]/count_fill[beam][evenodd][0];
          erlum_fill[beam][evenodd] = rlum_fill[beam][evenodd]*sqrt(e2count_fill[beam][evenodd][1]/pow2(count_fill[beam][evenodd][1])
              + e2count_fill[beam][evenodd][0]/pow2(count_fill[beam][evenodd][0]));
        }

      t_rlum->Fill();
    }


  f_rlum->cd();
  t_rlum->Write();
  t_rlum_check->Write();
  f_rlum->Close();
  f_bbc->Close();

  return 0;
}
