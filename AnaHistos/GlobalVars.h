/* pT bins */
const int npT = 30;
const double pTbin[npT+1] = { 0.0,
  0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
  5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
  12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

/* pT bins for ALL */
const int npT_pol = 15;
const double pTbin_pol[npT_pol+1] = { 2.0,
  2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0,
  10.0, 12.0, 15.0, 20.0, 30.0 };

/* Get ipt for TGraph Xaxis gx */
int Get_ipt(double *gx, double xpt)
{
  for(int ipt=0; ipt<npT; ipt++)
    if( TMath::Abs(gx[ipt] - xpt) < 0.2 )
      return ipt;

  cout << "Warning: No matching for pT = " << xpt << ", 0 returned!" << endl;
  return 0;
}

/* Get ipt_pol for given xpt */
int Get_ipt_pol(double xpt)
{
  for(int ipt=0; ipt<npT_pol; ipt++)
    if( xpt > pTbin_pol[ipt]-0.01 && xpt < pTbin_pol[ipt+1]+0.01 )
      return ipt;

  cout << "Warning: No matching for pT = " << xpt << ", 0 returned!" << endl;
  return 0;
}

TF1 *cross_pi0;
TF1 *cross_ph;
void SetWeight()
{
  // function for pT weight for pi0
  cross_pi0 = new TF1("cross_pi0", "x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))", 0.1, 100.);
  cross_pi0->SetParameters(2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00);

  // function for pT weight for direct photon
  cross_ph = new TF1("cross_ph", "x**(-[1]-[2]*log(x/[0]))*(1-(x/[0])**2)**[3]", 0.1, 100.);
  cross_ph->SetParameters(255., 5.98, 0.273, 14.43);
}
