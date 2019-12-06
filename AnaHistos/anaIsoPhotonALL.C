#include "GlobalVars.h"
#include "QueryTree.h"

void anaIsoPhotonALL(const int process = 0)
{
  const int nThread = 50;
  int thread = -1;

  TFile *f_rlum = new TFile("data/RelLum.root");
  TTree *t_rlum = (TTree*)f_rlum->Get("T");
  int runnumber, spin_pattern;
  double pol[2], epol[2], rlum[3][2], erlum[3][2];
  t_rlum->SetBranchAddress("Runnumber", &runnumber);
  t_rlum->SetBranchAddress("SpinPattern", &spin_pattern);
  t_rlum->SetBranchAddress("Pol", pol);
  t_rlum->SetBranchAddress("ePol", epol);
  t_rlum->SetBranchAddress("YellowRelLum", rlum[0]);
  t_rlum->SetBranchAddress("eYellowRelLum", erlum[0]);
  t_rlum->SetBranchAddress("BlueRelLum", rlum[1]);
  t_rlum->SetBranchAddress("eBlueRelLum", erlum[1]);
  t_rlum->SetBranchAddress("RelLum", rlum[2]);
  t_rlum->SetBranchAddress("eRelLum", erlum[2]);

  QueryTree *qt_ken2 = new QueryTree("data/YieldKEN2-isophoton.pdf");

  QueryTree *qt_asym = new QueryTree(Form("histos/isophoton-asym-%d.root",process), "RECREATE");

  for(int ien=0; ien<t_rlum->GetEntries(); ien++)
  {
    thread++;
    if( thread < process*nThread )
      continue;
    else if( thread >= (process+1)*nThread )
      break;

    t_rlum->GetEntry(ien);

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15717/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() )
    {
      cout << "Cannot open file for runnumber = " << runnumber << endl;
      break;
    }

    double nphoton[3][2][2][npT_pol] = {};  // beam, icr, ipol, ipt
    double npion[3][2][2][2][2][npT_pol] = {};  // beam, pttype(pt|2pt), ibg, icr, ipol, ipt

    int checkmap = 1;

    for(int beam=0; beam<3; beam++)
      for(int icr=0; icr<2; icr++)
        for(int ipol=0; ipol<2; ipol++)
        {
          int ih = beam + 3*icr + 3*2*ipol + 3*2*2*checkmap;
          TH1 *h_photon = (TH1*)f->Get(Form("h_1photon_pol_%d",ih));
          for(int ipt=0; ipt<npT_pol; ipt++)
          {
            int ptbin_first = h_photon->GetXaxis()->FindBin(pTbin_pol[ipt]);
            int ptbin_last = h_photon->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;
            nphoton[beam][icr][ipol][ipt] = h_photon->Integral(ptbin_first,ptbin_last);
          } // ipt

          for(int pttype=0; pttype<2; pttype++)
            for(int isoboth=0; isoboth<2; isoboth++)
              for(int isopair=0; isopair<2; isopair++)
              {
                char *ptname = pttype ? "2pt" : "";
                int ih = beam + 3*icr + 3*2*ipol + 3*2*2*checkmap + 3*2*2*2*isoboth + 3*2*2*2*2*isopair;
                TH2 *h2_pion = (TH2*)f->Get(Form("h2_2photon%s_pol_%d",ptname,ih));
                for(int ipt=0; ipt<npT_pol; ipt++)
                {
                  int ptbin_first = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt]);
                  int ptbin_last = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;
                  npion[beam][pttype][0][icr][ipol][ipt] += h2_pion->Integral(ptbin_first,ptbin_last, 113,162);
                  npion[beam][pttype][1][icr][ipol][ipt] += h2_pion->Integral(ptbin_first,ptbin_last, 48,97) +
                    h2_pion->Integral(ptbin_first,ptbin_last, 178,227);
                } // ipt
              } // pttype, isoboth, isopair
        } // beam, icr, ipol

    for(int beam=0; beam<3; beam++)
      for(int icr=0; icr<2; icr++)
      {
        double pbeam, e2pbeam;
        if(beam < 2)
        {
          pbeam = pol[beam];
          e2pbeam = epol[beam]*epol[beam]/pol[beam]/pol[beam];
        }
        else
        {
          pbeam = pol[0]*pol[1];
          e2pbeam = epol[0]*epol[0]/pol[0]/pol[0] + epol[1]*epol[1]/pol[1]/pol[1];
        }
        double r = rlum[beam][icr];
        double er = erlum[beam][icr];

        for(int ipt=0; ipt<npT_pol; ipt++)
        {

          for(int pttype=0; pttype<2; pttype++)
          {
            if( npion[beam][pttype][0][icr][1][ipt] + npion[beam][pttype][0][icr][0][ipt] < 10. ||
                npion[beam][pttype][0][icr][1][ipt] < 1. || npion[beam][pttype][0][icr][0][ipt] < 1. ||
                npion[beam][pttype][1][icr][1][ipt] < 1. || npion[beam][pttype][1][icr][0][ipt] < 1. )
              continue;

            for(int ibg=0; ibg<2; ibg++)
            {
              int imul = 1 + ibg + 2*pttype;
              int part = imul + 5*beam + 5*3*icr + 5*3*2*checkmap;
              double xpt, k2 = 1., ek2;
              //qt_ken2->Query(ipt, part, xpt, k2, ek2);

              double npp = npion[beam][pttype][ibg][icr][1][ipt];
              double npm = npion[beam][pttype][ibg][icr][0][ipt];

              double ALL = 1./pbeam*(npp - r*npm)/(npp + r*npm);
              double eALL = sqrt(pow(2*r*npp*npm/pbeam,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
                  + e2pbeam*ALL*ALL);

              if( TMath::Finite(ALL+eALL) && eALL > 0. )
              {
                int ig = imul + 5*beam + 5*3*icr + 5*3*2*spin_pattern + 5*3*2*4*ipt;
                qt_asym->Fill(runnumber, ig, runnumber, ALL, eALL);
              }
            } // ibg
          } // pttype

          if( nphoton[beam][icr][1][ipt] + nphoton[beam][icr][0][ipt] < 10. ||
              nphoton[beam][icr][1][ipt] < 1. || nphoton[beam][icr][0][ipt] < 1. )
            continue;

          int imul = 0;
          int part = imul + 5*beam + 5*3*icr + 5*3*2*checkmap;
          double xpt, k2 = 1., ek2;
          //qt_ken2->Query(ipt, part, xpt, k2, ek2);

          double npp = nphoton[beam][icr][1][ipt];
          double npm = nphoton[beam][icr][0][ipt];

          double ALL = 1./pbeam*(npp - r*npm)/(npp + r*npm);
          double eALL = sqrt(pow(2*r*npp*npm/pbeam,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
              + e2pbeam*ALL*ALL);

          if( TMath::Finite(ALL+eALL) && eALL > 0. )
          {
            int ig = imul + 5*beam + 5*3*icr + 5*3*2*spin_pattern + 5*3*2*4*ipt;
            qt_asym->Fill(runnumber, ig, runnumber, ALL, eALL);
          }
        } // ipt
      } // beam, icr

    delete f;
  } // ien

  qt_asym->Save();
}
