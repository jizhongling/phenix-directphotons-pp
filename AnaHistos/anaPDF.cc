// mysrc64 new
// g++ -std=c++11 -Wall -I$MYINSTALL/include -I$OFFLINE_MAIN/include -L$MYINSTALL/lib -L$OFFLINE_MAIN/lib `root-config --cflags --glibs` -lLHAPDF -o anaPDF anaPDF.cc
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "LHAPDF/LHAPDF.h"

using namespace std;

const int nrep = 1000;
const int npt = 7;
const int nx = 100;

const double pt[npt] = {6.423, 7.434, 8.443, 9.450, 10.83, 13.21, 16.89};
const double data[npt] = {0.002265, 0.002301, 0.004449, 0.003846, 0.01479, 0.01111, -0.02071};
const double err[npt] = {0.003100, 0.004800, 0.006900, 0.009500, 0.009901, 0.01400, 0.02200};

template<class T> inline constexpr T square(const T &x) { return x*x; }

void read_xsec(const char *fname, double xsec[][npt])
{
  ifstream fin(fname);
  char line[1024];

  while(fin.getline(line, 1024))
  {
    stringstream ss;
    string word;
    int irep, ipt;

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
      ss >> xsec_dir >> xsec_frag >> xsec[irep][ipt];
      ipt++;
    }
  }

  fin.close();
  return;
}

int main()
{
  double weight[nrep+1];
  double sumw = 0.;

  double unpol[1][npt];
  double pol[nrep+1][npt];
  read_xsec("data/dssv-unpol.txt", unpol);
  read_xsec("data/dssv-pol.txt", pol);

  for(int irep=1; irep<=nrep; irep++)
  {
    double chi2 = 0.;
    for(int ipt=0; ipt<npt; ipt++)
    {
      double all = pol[irep][ipt] / unpol[0][ipt];
      chi2 += square((all - data[ipt]) / err[ipt]);
    }
    weight[irep] = pow(chi2, (double)npt/2.-1.) * exp(-chi2/2.);
    sumw += weight[irep];
  }
  for(int irep=1; irep<=nrep; irep++)
    weight[irep] /= sumw;

  ofstream fout_old("data/reweighting-old.txt");
  ofstream fout_new("data/reweighting-new.txt");

  vector<LHAPDF::PDF*> v_pdf = LHAPDF::mkPDFs("DSSV_REP_LHAPDF6");

  for(int ix=0; ix<nx; ix++)
  {
    double x = 0.02 + 0.0006 * ix;
    double xg_center = 0.;
    double xg_old = 0.;
    double xg2_old = 0.;
    double xg_new = 0.;
    double xg2_new = 0.;

    for(int irep=0; irep<=nrep; irep++)
    {
      double xg = v_pdf.at(irep)->xfxQ2(21, x, 10.);
      if(irep == 0)
      {
        xg_center = xg;
      }
      else
      {
        xg_old += xg / (double)nrep;
        xg2_old += square(xg) / (double)nrep;
        xg_new += xg * weight[irep];
        xg2_new += square(xg) * weight[irep];
      }
    } // irep

    double exg_old = sqrt(xg2_old - square(xg_old));
    double exg_new = sqrt(xg2_new - square(xg_new));

    fout_old << x << "\t" << xg_center << "\t" << exg_old << endl;
    fout_new << x << "\t" << xg_new << "\t" << exg_new << endl;
  }

  fout_old.close();
  fout_new.close();

  return 0;
}
