// mysrc64 new
// g++ -std=c++11 -Wall -I$MYINSTALL/include -L$MYINSTALL/lib -lLHAPDF -o anaPDF_dgsgn anaPDF_dgsgn.cc
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "LHAPDF/LHAPDF.h"

using namespace std;

const int npt = 7 + 2;
const double pt[npt] = {5., 6.423, 7.434, 8.443, 9.450, 10.83, 13.21, 16.89, 20.};
const double data[npt] = {0., -0.002265, 0.002301, 0.004449, 0.003846, 0.01479, 0.01111, -0.02071, 0.};
const double err[npt] = {0., 0.003100, 0.004800, 0.006900, 0.009500, 0.009901, 0.01400, 0.02200, 0.};

template<class T> inline constexpr T square(const T &x) { return x*x; }

void read_xsec(const char *fname, double xsec[][npt], const int iadd = 0)
{
  ifstream fin(fname);
  char line[1024];
  int irep, ipt;

  while(fin.getline(line, 1024))
  {
    stringstream ss;
    string word;

    ss << line;
    ss >> word;
    if(word.compare(0, 7, "REPLICA") == 0)
    {
      ss >> irep;
      ipt = 0;
    }
    else if(!word.empty())
    {
      double xsec_dir, xsec_frag;
      ss >> xsec_dir >> xsec_frag >> xsec[irep+iadd][ipt];
      ipt++;
    }
  }

  fin.close();
  return;
}

int main()
{
  const int nrep_pos = 461;
  const int nrep_neg = 72;
  const int nrep = nrep_pos + nrep_neg;

  double unpol[1][npt];
  double pol[nrep+1][npt];
  read_xsec("data/cross-unpol-NNPDF30_nlo_as_0118.txt", unpol);
  read_xsec("data/cross-pol-JAM22_pol_SU23_pos_g.txt", pol);
  read_xsec("data/cross-pol-JAM22_pol_SU23_neg_g.txt", pol, nrep_pos);

  int irep_min = 9999;
  double chi2_min = 9999.;
  for(int irep=1; irep<=nrep; irep++)
  {
    double chi2 = 0.;
    for(int ipt=1; ipt<npt-1; ipt++)
    {
      double all = pol[irep][ipt] / unpol[0][ipt];
      chi2 += square((all - data[ipt]) / err[ipt]);
    }
    if(chi2 < chi2_min)
    {
      irep_min = irep;
      chi2_min = chi2;
    }
  }

  double weight[nrep+1];
  for(int irep=1; irep<=nrep; irep++)
  {
    double chi2 = 0.;
    for(int ipt=1; ipt<npt-1; ipt++)
    {
      double all = pol[irep][ipt] / unpol[0][ipt];
      chi2 += square((all - data[ipt]) / err[ipt]);
    }
    // See Erratum of Nucl. Phys. B 849 (2011) 112-143
    weight[irep] = pow(chi2, (npt-2-1)/2.) * exp(-chi2/2.);
  }

  ofstream fout_dir("data/all-JAM22_pol_SU23-chi2min.txt");
  ofstream fout_jam[3];
  fout_jam[0].open("data/all-JAM22_pol_SU23-dgpos.txt");
  fout_jam[1].open("data/all-JAM22_pol_SU23-dgneg.txt");
  fout_jam[2].open("data/all-JAM22_pol_SU23-reweight.txt");
  //vector<LHAPDF::PDF*> v_pdf = LHAPDF::mkPDFs("JAM22ppdf");

  for(int ipt=0; ipt<npt; ipt++)
  {
    double sum_all[3] = {};
    double sum_all2[3] = {};
    double n_all[3] = {};

    for(int irep=1; irep<=nrep; irep++)
    {
      //int ixg = v_pdf.at(irep)->xfxQ2(21, 0.5, 10.) < 0 ? 1 : 0;
      int ixg = irep > nrep_pos ? 1 : 0;
      double all = pol[irep][ipt] / unpol[0][ipt];
      sum_all[ixg] += all;
      sum_all2[ixg] += square(all);
      n_all[ixg]++;
      sum_all[2] += all * weight[irep];
      sum_all2[2] += square(all) * weight[irep];
      n_all[2] += weight[irep];
      if(irep == irep_min)
        fout_dir << pt[ipt] << "\t" << all << endl;
    }

    for(int ixg=0; ixg<3; ixg++)
    {
      double mean_all = sum_all[ixg] / n_all[ixg];
      double sigma_all = sqrt(sum_all2[ixg] / n_all[ixg] - square(mean_all));
      fout_jam[ixg] << pt[ipt] << "\t" << mean_all << "\t" << sigma_all << endl;
    }
  }

  fout_dir.close();
  fout_jam[0].close();
  fout_jam[1].close();
  fout_jam[2].close();

  return 0;
}
