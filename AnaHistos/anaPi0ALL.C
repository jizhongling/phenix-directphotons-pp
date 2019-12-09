#include "GlobalVars.h"
#include "QueryTree.h"

void anaPi0ALL()
{
  TFile *f_rlum = new TFile("data/RelLum.root");
  TTree *t_rlum = (TTree*)f_rlum->Get("T");
  int runnumber, spin_pattern;
  double pol[2], epol[2], rlum[2], erlum[2];
  t_rlum->SetBranchAddress("Runnumber", &runnumber);
  t_rlum->SetBranchAddress("SpinPattern", &spin_pattern);
  t_rlum->SetBranchAddress("Pol", pol);
  t_rlum->SetBranchAddress("ePol", epol);
  t_rlum->SetBranchAddress("RelLum", rlum);
  t_rlum->SetBranchAddress("eRelLum", erlum);

  QueryTree *qt_ken2 = new QueryTree("data/YieldKEN2-pion.pdf");

  TFile *f_ken2 = new TFile("/phenix/plhf/zji/sources/offline/AnalysisTrain/Run13_Pi0Ana_YIS/KENFactor.root");
  TTree *t_ken2 = (TTree*)f_ken2->Get("T");
  double k_total_even, k_total_odd, k_back_even, k_back_odd;
  t_ken2->SetBranchAddress("K_Total_Even", &k_total_even);
  t_ken2->SetBranchAddress("K_Total_Odd", &k_total_odd);
  t_ken2->SetBranchAddress("K_Back_Even", &k_back_even);
  t_ken2->SetBranchAddress("K_Back_Odd", &k_back_odd);

  double ken2[2][2][npT_pol];  // ibg, icr, ipt
  for(int ibg=0; ibg<2; ibg++)
    for(int icr=0; icr<2; icr++)
      for(int ipt=0; ipt<npT_pol; ipt++)
        ken2[ibg][icr][ipt] = 1.;

  for(int ien=0; ien<t_ken2->GetEntries(); ien++)
  {
    t_ken2->GetEntry(ien);
    ken2[0][0][ien] = k_total_even;
    ken2[0][1][ien] = k_total_odd;
    ken2[1][0][ien] = k_back_even;
    ken2[1][1][ien] = k_back_odd;
  }

  QueryTree *qt_asym = new QueryTree("data/pion-asym.root", "RECREATE");

  for(int ien=0; ien<t_rlum->GetEntries(); ien++)
  {
    if(ien && ien%100 == 0)
      cout << ien << endl;

    t_rlum->GetEntry(ien);

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15673/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() )
    {
      cout << "Cannot open file for runnumber = " << runnumber << endl;
      continue;
    }

    double npion[2][2][2][npT_pol];  // ibg, icr, ipol, ipt

    for(int icr=0; icr<2; icr++)
      for(int ipol=0; ipol<2; ipol++)
      {
        int ih = icr + 2*ipol;
        TH2 *h2_pion = (TH2*)f->Get(Form("h2_pion_pol_%d",ih));

        for(int ipt=0; ipt<npT_pol; ipt++)
        {
          int ptbin_first = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt]);
          int ptbin_last = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;

          npion[0][icr][ipol][ipt] = h2_pion->Integral(ptbin_first,ptbin_last, 113,162);
          npion[1][icr][ipol][ipt] = h2_pion->Integral(ptbin_first,ptbin_last, 48,97) +
            h2_pion->Integral(ptbin_first,ptbin_last, 178,227);
        } // ipt
      } // icr, ipol

    for(int ibg=0; ibg<2; ibg++)
      for(int icr=0; icr<2; icr++)
        for(int ipt=0; ipt<npT_pol; ipt++)
        {
          if( npion[0][icr][1][ipt] + npion[0][icr][0][ipt] < 10. ||
              npion[0][icr][1][ipt] < 1. || npion[0][icr][0][ipt] < 1. ||
              npion[1][icr][1][ipt] < 1. || npion[1][icr][0][ipt] < 1. )
            continue;

          double pb = pol[0];
          double py = pol[1];
          double epb = epol[0];
          double epy = epol[1];
          double r = rlum[icr];
          double er = erlum[icr];

          int part = ibg + 2*icr;
          double xpt, k2, ek2;
          qt_ken2->Query(ipt, part, xpt, k2, ek2);
          //k2 = ken2[ibg][icr][ipt];

          double npp = npion[ibg][icr][1][ipt];
          double npm = npion[ibg][icr][0][ipt];

          double ALL = 1./(pb*py)*(npp - r*npm)/(npp + r*npm);
          double eALL = sqrt(pow(2*r*npp*npm/pb/py,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
              + (epb*epb/pb/pb + epy*epy/py/py)*ALL*ALL);

          if( TMath::Finite(ALL+eALL) && eALL > 0. )
          {
            int ig = ibg + 2*icr + 2*2*spin_pattern + 2*2*4*ipt;
            qt_asym->Fill(runnumber, ig, runnumber, ALL, eALL);
          }
        } // ibg, icr, ipt

    delete f;
  } // ien

  qt_asym->Save();
}
