/*
 * pT bins
 */
const Int_t npT = 30;
const Double_t pTbin[npT+1] = { 0.0,
  0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
  5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
  12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

/* Get ipt for TGraph Xaxis gx */
Int_t Get_ipt(Double_t *gx, Double_t xx)
{
  for(Int_t ipt=0; ipt<npT; ipt++)
    if( TMath::Abs(gx[ipt] - xx) < 0.2 )
      return ipt;

  cout << "Warning: No matching for pT = " << xx << ", 0 returned!" << endl;
  return 0;
}
