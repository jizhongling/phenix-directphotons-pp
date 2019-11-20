#include "GlobalVars.h"
#include "QueryTree.h"

void anaIsoPhotonALL()
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

  QueryTree *qt_ken2 = new QueryTree("data/YieldKEN2-isophoton.pdf");

  QueryTree *qt_asym = new QueryTree("data/isophoton-asym.root", "RECREATE");

  for(int ien=0; ien<t_rlum->GetEntries(); ien++)
  {
    if(ien && ien%100 == 0)
      cout << ien << endl;

    t_rlum->GetEntry(ien);

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15666/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() )
    {
      cout << "Cannot open file for runnumber = " << runnumber << endl;
      break;
    }

    double nphoton[2][2][npT_pol] = {};  // icr, ipol, ipt
    double npion[2][2][2][2][2][npT_pol] = {};  // pttype(pt|2pt), isotype(isoboth|isopair), ibg, icr, ipol, ipt

    int checkmap = 1;

    for(int icr=0; icr<2; icr++)
      for(int ipol=0; ipol<2; ipol++)
      {
        int ih = icr + 2*ipol + 2*2*checkmap;
        TH1 *h_photon = (TH1*)f->Get(Form("h_1photon_pol_%d",ih));
        for(int ipt=0; ipt<npT_pol; ipt++)
        {
          int ptbin_first = h_photon->GetXaxis()->FindBin(pTbin_pol[ipt]);
          int ptbin_last = h_photon->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;
          nphoton[icr][ipol][ipt] = h_photon->Integral(ptbin_first,ptbin_last);
        } // ipt

        for(int ptype=0; pttype<2; pttype++)
          for(int isotype=0; isotype<2; isotype++)
            for(int iso=0; iso<2; iso++)
            {
              const char *ptname = pttype ? "2pt" : "";
              int ih = icr + 2*ipol + 2*2*checkmap + 2*2*2*(1-isotype*iso) + 2*2*2*2*(1-(1-isotype)*iso);
              TH2 *h2_pion = (TH2*)f->Get(Form("h2_2photon%s_pol_%d",ptname,ih));
              for(int ipt=0; ipt<npT_pol; ipt++)
              {
                int ptbin_first = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt]);
                int ptbin_last = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;
                npion[pttype][isotype][0][icr][ipol][ipt] += h2_pion->Integral(ptbin_first,ptbin_last, 113,162);
                npion[pttype][isotype][1][icr][ipol][ipt] += h2_pion->Integral(ptbin_first,ptbin_last, 48,97) +
                  h2_pion->Integral(ptbin_first,ptbin_last, 178,227);
              } // ipt
            } // pttype, isotype, iso
      } // icr, ipol

    for(int icr=0; icr<2; icr++)
      for(int ipt=0; ipt<npT_pol; ipt++)
      {
        double pb = pol[0];
        double py = pol[1];
        double epb = epol[0];
        double epy = epol[1];
        double r = rlum[icr];
        double er = erlum[icr];

        for(int ptype=0; pttype<2; pttype++)
          for(int isoype=0; isotype<2; isotype++)
          {
            if( npion[pttype][isotype][0][icr][1][ipt] + npion[pttype][isotype][0][icr][0][ipt] < 10. ||
                npion[pttype][isotype][0][icr][1][ipt] < 1. || npion[pttype][isotype][0][icr][0][ipt] < 1. ||
                npion[pttype][isotype][1][icr][1][ipt] < 1. || npion[pttype][isotype][1][icr][0][ipt] < 1. )
              continue;

            for(int ibg=0; ibg<2; ibg++)
            {
              int part = pttype + 2*isotype + 2*2*ibg + 2*2*2*icr;
              double xpt, k2, ek2;
              qt_ken2->Query(ipt, part, xpt, k2, ek2);

              double npp = npion[pttype][isotype][ibg][icr][1][ipt];
              double npm = npion[pttype][isotype][ibg][icr][0][ipt];

              double ALL = 1./(pb*py)*(npp - r*npm)/(npp + r*npm);
              double eALL = sqrt(pow(2*r*npp*npm/pb/py,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
                  + (epb*epb/pb/pb + epy*epy/py/py)*ALL*ALL);

              if( TMath::Finite(ALL+eALL) && eALL > 0. )
              {
                int ig = ibg + 3*icr + 3*2*spin_pattern + 3*2*4*ipt;
                qt_asym->Fill(runnumber, ig, runnumber, ALL, eALL);
              }
            } // ibg
          } // pttype, isotype

        if( nphoton[icr][1][ipt] + nphoton[icr][0][ipt] < 10. ||
            nphoton[icr][1][ipt] < 1. || nphoton[icr][0][ipt] < 1. )
          continue;

        int part = 2*2*2*2 + icr;
        double xpt, k2, ek2;
        qt_ken2->Query(ipt, part, xpt, k2, ek2);

        double npp = nphoton[icr][1][ipt];
        double npm = nphoton[icr][0][ipt];

        double ALL = 1./(pb*py)*(npp - r*npm)/(npp + r*npm);
        double eALL = sqrt(pow(2*r*npp*npm/pb/py,2)/pow(npp + r*npm,4)*(k2/npp + k2/npm + er*er/r/r)
            + (epb*epb/pb/pb + epy*epy/py/py)*ALL*ALL);

        if( TMath::Finite(ALL+eALL) && eALL > 0. )
        {
          int ig = 2 + 3*icr + 3*2*spin_pattern + 3*2*4*ipt;
          qt_asym->Fill(runnumber, ig, runnumber, ALL, eALL);
        }
      } // icr, ipt

    delete f;
  } // ien

  qt_asym->Save();
}
