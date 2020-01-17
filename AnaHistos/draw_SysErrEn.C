#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"

void draw_SysErrEn()
{
  const double corr[4] = {1., 0.993, 1., 1.};
  const double A = 0.35;

  QueryTree *qt_sys = new QueryTree("data/syserr-en.root", "RECREATE");

  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-syserr.root");
  THnSparse *hn_photon = (THnSparse*)f->Get("hn_photon");
  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  hn_photon->GetAxis(4)->SetRange(3,3); // econe_trk[ival]: EMCal, nomap, withmap
  hn_pion->GetAxis(6)->SetRange(3,3); // econe_trk[ival]: EMCal, nomap, withmap

  for(int ipt=10; ipt<npT; ipt++)
  {
    double ndir[2][4], endir[2][4];  // isolated, isys
    double xpt, Miss, eMiss;
    qt_miss->Query(ipt, 3, xpt, Miss, eMiss);

    for(int isys=0; isys<4; isys++)
    {
      hn_photon->GetAxis(5)->SetRange(isys+1,isys+1);
      hn_pion->GetAxis(7)->SetRange(isys+1,isys+1);

      double nphoton[2], enphoton[2];  // isolated
      for(int isolated=0; isolated<2; isolated++)
      {
        hn_photon->GetAxis(3)->SetRange(isolated+1,2);
        TH1 *h_1photon = hn_photon->Projection(0);  // pt_truth
        nphoton[isolated] = h_1photon->GetBinContent(ipt+1);
        enphoton[isolated] = h_1photon->GetBinError(ipt+1);
        delete h_1photon;
      }

      double npion[2][2], enpion[2][2];  // isoboth, isopair
      for(int isoboth=0; isoboth<2; isoboth++)
        for(int isopair=0; isopair<2; isopair++)
          if(isoboth != 1 || isopair != 1)
          {
            TH1 *h_minv;
            hn_pion->GetAxis(4)->SetRange(isoboth+1,2);  // isoboth
            hn_pion->GetAxis(5)->SetRange(isopair+1,2);  // isopair
            hn_pion->GetAxis(0)->SetRange(ipt+1,ipt+1);  // pt_truth
            h_minv = (TH1*)hn_pion->Projection(2);  // minv
            FitMinv(h_minv, npion[isoboth][isopair], enpion[isoboth][isopair], false, 0.11*corr[isys],0.16*corr[isys]);
            delete h_minv;
          }

      ndir[0][isys] = nphoton[0] - (1 + Miss)*(1 + A)*npion[0][0];
      endir[0][isys] = sqrt(pow(enphoton[0],2) + pow((1 + Miss)*(1 + A)*enpion[0][0],2));

      ndir[1][isys] = nphoton[1] - npion[1][0] - Miss*(1 + A)*npion[0][1];
      endir[1][isys] = sqrt(pow(enphoton[1],2) + pow(enpion[1][0],2) + pow(Miss*(1 + A)*enpion[0][1],2));

      if(isys > 0)
        for(int isolated=0; isolated<2; isolated++)
        {
          int index = isolated + 2*(isys-1);
          double rsys = ndir[isolated][isys]/ndir[isolated][0];
          double ersys = rsys*sqrt(pow(endir[isolated][isys]/ndir[isolated][isys],2) + pow(endir[isolated][0]/ndir[isolated][0],2));
          if( TMath::Finite(rsys+ersys) )
            qt_sys->Fill(ipt, index, xpt, rsys, ersys);
        } // isolated
    } // isys
  } // ipt

  qt_sys->Save();
}
