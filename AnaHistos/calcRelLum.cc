// To compile: g++ -Wall -o calcRelLum calcRelLum.cc -I$OFFLINE_MAIN/include -L$OFFLINE_MAIN/lib -m32 -luspin -lodbc -lgsl -lgslcblas -lm `root-config --cflags --libs`
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/foreach.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <TFile.h>
#include <TTree.h>

#include <SpinDBOutput.hh>
#include <SpinDBContent.hh>

using namespace std;

const double kn = 0.2278;
const double ks = 0.2278;

double rate_f(double x, void *params)
{
  double *p = (double*)params;
  double rate_obs = *p;

  double rate_true = 1 - exp(-x*(1+kn)) - exp(-x*(1+ks)) + exp(-x*(1+kn+ks));

  return rate_true - rate_obs;
}

double rate_df(double x, void *params)
{
  double df = (1+kn)*exp(-x*(1+kn)) + (1+ks)*exp(-x*(1+ks)) - (1+kn+ks)*exp(-x*(1+kn+ks));

  return df;
}

void rate_fdf(double x, void *params, double *y, double *dy)
{
  double *p = (double*)params;
  double rate_obs = *p;

  double expkn = exp(-x*(1+kn));
  double expks = exp(-x*(1+ks));
  double expknks = exp(-x*(1+kn+ks));

  double rate_true = 1 - expkn - expks + expknks;
  double df = (1+kn)*expkn + (1+ks)*expks - (1+kn+ks)*expknks;

  *y = rate_true - rate_obs;
  *dy = df;
}

// Input: rate_obs[2] = {rate_vtx, rate_novtx}
// Output: rate_true[3] = {rate_vtx_true, rate_novtx_true, rate_vtx_res}
void rate_corr(double rate_obs[], double rate_true[])
{
  const int max_iter = 100;
  const double epsrel = 1e-6;

  const gsl_root_fdfsolver_type *solver_type = gsl_root_fdfsolver_newton;
  gsl_root_fdfsolver *solver = gsl_root_fdfsolver_alloc(solver_type);

  double params;
  gsl_function_fdf FDF;
  FDF.f = &rate_f;
  FDF.df = &rate_df;
  FDF.fdf = &rate_fdf;
  FDF.params = &params;

  for(int i=0; i<3; i++)
  {
    if(i < 2)
      params = rate_obs[i];
    else
      params = rate_obs[0]*(rate_obs[0]/rate_obs[1])/(rate_true[0]/rate_true[1]);

    int status;
    int iter = 0;
    double x0, x = params;
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
  }

  gsl_root_fdfsolver_free(solver);
  return;
}

void get_rlum_gl1p(SpinDBOutput &spin_out, SpinDBContent &spin_cont,
    int runnumber, double pol[], double epol[], double rlum[])
{
  /* Retrieve entry from Spin DB and get fill number */
  int qa_level = spin_out.GetDefaultQA(runnumber);
  spin_out.StoreDBContent(runnumber, runnumber, qa_level);
  spin_out.GetDBContentStore(spin_cont, runnumber);

  if( spin_out.CheckRunRow(runnumber,qa_level) != 1 )
  {
    cerr << "Wrong spinDB info!" << endl;
    exit(1);
  }

  double pbsyst;
  double pysyst;
  spin_cont.GetPolarizationBlue(1, pol[0], epol[0], pbsyst);
  spin_cont.GetPolarizationYellow(1, pol[1], epol[1], pysyst);

  long long bbc30_gl1p[2][2] = {};

  for(int ib=0; ib<120; ib++)
  {
    int pb = spin_cont.GetSpinPatternBlue(ib);
    int py = spin_cont.GetSpinPatternYellow(ib);
    int pattern = (pb*py + 1)/2;
    if(pattern == 0 || pattern == 1)
      bbc30_gl1p[ib%2][pattern] += spin_cont.GetScaler(0, ib);
  }

  for(int evenodd=0; evenodd<2; evenodd++)
    rlum[evenodd] = (double)bbc30_gl1p[evenodd][1]/bbc30_gl1p[evenodd][0];

  return;
}

int main()
{
  int runnumber;
  int spin_pattern = 0;

  TFile *f_rlum = new TFile("data/RelLum.root", "RECREATE");

  TTree *t_rlum = new TTree("T", "Relative luminosity");
  double pol[2], epol[2], rlum[2];
  t_rlum->Branch("Runnumber", &runnumber, "Runnumber/I");
  t_rlum->Branch("SpinPattern", &spin_pattern, "SpinPattern/I");
  t_rlum->Branch("Pol", pol, "Pol[2]/D");
  t_rlum->Branch("ePol", epol, "ePol[2]/D");
  t_rlum->Branch("RelLum", rlum, "RelLum[2]/D");

  TTree *t_rlum_check = new TTree("T1", "Relative luminosity check");
  double gl1p[2], ss_raw[2], ss_pile[2], ss_res[2];
  t_rlum_check->Branch("Runnumber", &runnumber, "Runnumber/I");
  t_rlum_check->Branch("GL1p", gl1p, "Gl1p[2]/D");
  t_rlum_check->Branch("SS_Uncorr", ss_raw, "SS_Uncorr[2]/D");
  t_rlum_check->Branch("SS_Pileup", ss_pile, "SS_Pileup[2]/D");
  t_rlum_check->Branch("SS_Residual", ss_res, "SS_Residual[2]/D");

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
  int pb[120], py[120], clock[120], bbcn[120], bbcs[120],
      bbcns[120], bbc30[120], bbcnovtx[120];
  t_bbc->SetBranchAddress("Runnumber", &runnumber);
  t_bbc->SetBranchAddress("SpinPatBlue", pb);
  t_bbc->SetBranchAddress("SpinPatYellow", py);
  t_bbc->SetBranchAddress("CLOCK", clock);
  t_bbc->SetBranchAddress("BBCN", bbcn);
  t_bbc->SetBranchAddress("BBCS", bbcs);
  t_bbc->SetBranchAddress("BBCNS", bbcns);
  t_bbc->SetBranchAddress("BBC_30cm", bbc30);
  t_bbc->SetBranchAddress("BBCnovtx", bbcnovtx);

  vector<int> runnoDone(1000);

  for(int ien=0; ien<t_bbc->GetEntries(); ien++)
  {
    t_bbc->GetEntry(ien);

    get_rlum_gl1p(spin_out, spin_cont, runnumber, pol, epol, gl1p);

    double rate_raw[2][2] = {};
    double rate_pile[2][2] = {};
    double rate_lum[2][2] = {};

    for(int ib=0; ib<120; ib++)
      if( abs(pb[ib]*py[ib]) == 1 && clock[ib] > 0 && bbcn[ib] > 0 && bbcs[ib] > 0 &&
          bbcns[ib] > 0 && bbc30[ib] > 0 && bbcnovtx[ib] > 0 )
      {
        int pattern = (pb[ib]*py[ib] + 1)/2;

        double bbc_vtx = (double)bbc30[ib]/clock[ib];
        double bbc_novtx = (double)bbcnovtx[ib]/clock[ib];
        double rate_obs[2] = {bbc_vtx, bbc_novtx};
        double rate_true[3];
        rate_corr(rate_obs, rate_true);

        rate_raw[ib%2][pattern] += rate_obs[0];
        rate_pile[ib%2][pattern] += rate_true[0];
        rate_lum[ib%2][pattern] += rate_true[2];
      }

    for(int evenodd=0; evenodd<2; evenodd++)
    {
      ss_raw[evenodd] = rate_raw[evenodd][1]/rate_raw[evenodd][0];
      ss_pile[evenodd] = rate_pile[evenodd][1]/rate_pile[evenodd][0];
      ss_res[evenodd] = rate_lum[evenodd][1]/rate_lum[evenodd][0];
      rlum[evenodd] = ss_res[evenodd];
    }

    map_int_t::iterator it_Inseok = runnoInseok.find(runnumber);
    if( it_Inseok != runnoInseok.end() )
    {
      spin_pattern = it_Inseok->second;

      if(runnumber < 386946)
        for(int evenodd=0; evenodd<2; evenodd++)
          rlum[evenodd] = gl1p[evenodd];

      t_rlum->Fill();
      runnoDone.push_back(runnumber);
    }

    t_rlum_check->Fill();
  }

  BOOST_FOREACH(const map_int_t::value_type &runno, runnoInseok)
    if( find(runnoDone.begin(), runnoDone.end(), runno.first) == runnoDone.end() )
    {
      runnumber = runno.first;
      spin_pattern = runno.second;

      get_rlum_gl1p(spin_out, spin_cont, runnumber, pol, epol, rlum);

      t_rlum->Fill();
    }


  f_rlum->cd();
  t_rlum->Write();
  t_rlum_check->Write();
  f_rlum->Close();
  f_bbc->Close();

  return 0;
}
