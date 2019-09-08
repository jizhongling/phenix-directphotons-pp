#include "PtWeights.h"

#include <TOAD.h>

#include <TString.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2.h>

using namespace std;

PtWeights::PtWeights():
  cross_pi0(nullptr),
  cross_ph(nullptr),
  f_mb(nullptr),
  h2_mb(nullptr)
{
  /* Function for pT weights for pi0 */
  cross_pi0 = new TF1("cross_pi0", "x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))", 0.1, 100.);
  if(!cross_pi0)
  {
    cerr << "No cross_pi0" << endl;
    exit(1);
  }
  cross_pi0->SetParameters(2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00);

  /* Function for pT weights for direct photon */
  cross_ph = new TF1("cross_ph", "x**(-[1]-[2]*log(x/[0]))*(1-(x/[0])**2)**[3]", 0.1, 100.);
  if(!cross_ph)
  {
    cerr << "No cross_ph" << endl;
    exit(1);
  }
  cross_ph->SetParameters(255., 5.98, 0.273, 14.43);

  ReadWeights();
}

PtWeights::~PtWeights()
{
  delete cross_pi0;
  delete cross_ph;
  delete f_mb;
}

double PtWeights::EvalPi0(double pt)
{
  return cross_pi0->Eval(pt);
}

double PtWeights::EvalPhoton(double pt)
{
  return cross_ph->Eval(pt);
}

double PtWeights::Integral(double pt1, double pt2, const char *option)
{
  double weight = 1.;

  TString opt = option;
  opt.ToLower();
  if( opt.Contains("pion") )
  {
    weight = cross_pi0->Integral(pt1, pt2);
  }
  else if( opt.Contains("photon") )
  {
    weight = cross_ph->Integral(pt1, pt2);
  }
  else if( opt.Contains("minbias") )
  {
    int binx1 = h2_mb->GetXaxis()->FindBin(pt1);
    int binx2 = h2_mb->GetXaxis()->FindBin(pt2);

    weight = h2_mb->Integral(binx1, binx2);
  }

  return weight;
}

void PtWeights::ReadWeights()
{
  TOAD *toad_loader = new TOAD("AnaFastMC");
  toad_loader->SetVerbosity(0);
  string file_location = toad_loader->location("MinBiasPtWeights.root");
  cout << "TOAD file location: " << file_location << endl;

  f_mb = new TFile( file_location.c_str() );
  if(!f_mb || f_mb->IsZombie())
  {
    cout << "Cannot open " << file_location << endl;
    exit(1);
  }
  h2_mb = (TH2*)f_mb->Get("h2_proc_pt");
  if(!h2_mb)
  {
    cout << "No h2_mb" << endl;
    exit(1);
  }

  delete toad_loader;
}
