const int ngr_pi0 = 3*2*4*2;
const int ngr_pion = 3*2*2*4*2;
const int ngr_photon = 3*2*4;

const double pt_eta[6] = {2.46, 3.35, 4.38, 5.39, 6.41, 7.69};
const double all_eta[6] = {-27, 1, 2, 100, 80, 130};
const double stat_eta[6] = {37, 44, 77, 140, 230, 410};
const double syst_eta[6] = {4.6, 4.7, 4.8, 5.5, 5.0, 11};

double Get_sys_eta(double xpt)
{
  const double scale = 510./200.;
  int iptmin = 0;
  double diffmin = 9999.;
  for(int ipt=0; ipt<6; ipt++)
  {
    double diff = fabs(xpt - pt_eta[ipt]*scale);
    if(diff < diffmin)
    {
      iptmin = ipt;
      diffmin = diff;
    }
  }

  double sys = 1e-4*sqrt(stat_eta[iptmin]*stat_eta[iptmin] + syst_eta[iptmin]*syst_eta[iptmin]);
  return sys;
}
