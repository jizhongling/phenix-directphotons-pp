const int npT = 7;
const int pTlow[2][npT] = { {13, 13, 14, 15, 16, 17, 18}, {5, 3, 4, 5, 6, 7, 8} };
const int pThigh[2][npT] = { {30, 14, 15, 16, 17, 18, 30}, {20, 4, 5, 6, 7, 8, 30} };
const double pTbin[2][npT] = { {6., 6.5, 7., 7.5, 8., 8.5, 9.5}, {1., 1.5, 2., 2.5, 3., 3.5, 4.5} };
const double pTbinL[2][npT] = { {6., 6., 6.5, 7., 7.5, 8., 8.5}, {2., 1., 1.5, 2., 2.5, 3., 3.5} };
const double pTbinC[2][npT] = { {6.25, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75}, {2.25, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75} };

/* Get ipt for TGraph Xaxis gx */
int Get_ipt(double *gx, double xx)
{
  for(int ipt=0; ipt<30; ipt++)
    if( TMath::Abs(gx[ipt] - xx) < 0.2 )
      return ipt;

  cout << "Warning: No matching for pT = " << xx << ", 0 returned!" << endl;
  return 0;
}

