#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_SysErrEn()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  const double Prob = 0.96;
  const double eProb = 0.02;
  const double ToF[3] = {0.992, 0.992, 0.997};
  const double eToF[3] = {0.002, 0.002, 0.002};
  const double Conv[3] = {0.849, 0.971, 0.959};
  const double eConv[3] = {0.027, 0.023, 0.029};
  const double A = 0.28;
  const double eA = 0.05;

  QueryTree *qt_sys = new QueryTree("data/syserr-en.root", "RECREATE");

  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");
  QueryTree *qt_miss_eta = new QueryTree("data/MissingRatio-eta.root");
  QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root");
  QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root");
  QueryTree *qt_badpass = new QueryTree("data/MergePassRate.root");
  QueryTree *qt_veto = new QueryTree("data/SelfVeto.root");
  qt_merge1->SetQuiet();
  qt_merge2->SetQuiet();

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-syserr.root");
  THnSparse *hn_photon = (THnSparse*)f->Get("hn_photon");
  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  hn_photon->GetAxis(4)->SetRange(3,3); // econe_trk[ival]: EMCal, nomap, withmap
  hn_pion->GetAxis(6)->SetRange(3,3); // econe_trk[ival]: EMCal, nomap, withmap

  for(int isoboth=0; isoboth<2; isoboth++)
    for(int isopair=0; isopair<2; isopair++)
      for(int pttype=0; pttype<2; pttype++)
        if( (1-isoboth)*(1-isopair) == 1 ||
            isoboth*(1-isopair)*(1-pttype) == 1 ||
            (1-isoboth)*isopair == 1 )
          for(int part=0; part<3; part++)
            for(int isys=0; isys<4; isys++)
            {
              int ic = isoboth + 2*isopair + 2*2*pttype + 2*2*2*part + 2*2*2*3*isys;
              mc(ic, 6,5);
            }

  for(int ipt=10; ipt<npT; ipt++)
  {
    double xpt, nsys[2][4], ensys[2][4];  // isolated, isys

    for(int isys=0; isys<4; isys++)
    {
      double ndir[2][3], endir[2][3];  // isolated, part
      hn_photon->GetAxis(5)->SetRange(isys+1,isys+1);
      hn_pion->GetAxis(7)->SetRange(isys+1,isys+1);

      for(int part=0; part<3; part++)
      {
        double nph[2], enph[2];  // isolated
        for(int isolated=0; isolated<2; isolated++)
        {
          hn_photon->GetAxis(3)->SetRange(isolated+1,2);
          hn_photon->GetAxis(2)->SetRange(secl[part],sech[part]);  // sector
          TH1 *h_1photon = hn_photon->Projection(1);  // ptsim
          nph[isolated] = h_1photon->GetBinContent(ipt+1);
          enph[isolated] = h_1photon->GetBinError(ipt+1);
          delete h_1photon;
        }

        double npion[2][2][2], enpion[2][2][2];  // isoboth, isopair, pttype
        for(int isoboth=0; isoboth<2; isoboth++)
          for(int isopair=0; isopair<2; isopair++)
            for(int pttype=0; pttype<2; pttype++)
              if( (1-isoboth)*(1-isopair) == 1 ||
                  isoboth*(1-isopair)*(1-pttype) == 1 ||
                  (1-isoboth)*isopair == 1 )
              {
                int ic = isoboth + 2*isopair + 2*2*pttype + 2*2*2*part + 2*2*2*3*isys;
                mcd(ic, ipt+1);
                TH1 *h_minv;
                hn_pion->GetAxis(4)->SetRange(isoboth+1,2);  // isoboth
                hn_pion->GetAxis(5)->SetRange(isopair+1,2);  // isopair
                hn_pion->GetAxis(3)->SetRange(secl[part],sech[part]);  // sector
                hn_pion->GetAxis(pttype)->SetRange(ipt+1,ipt+1);  // 0: pt_1photon
                hn_pion->GetAxis(1-pttype)->SetRange(1,npT);  // 1: pt_2photon
                h_minv = (TH1*)hn_pion->Projection(2);  // minv
                h_minv->Rebin(10);
                h_minv->SetTitle( Form("part %d p_{T}: %3.1f-%3.1f GeV", part, pTbin[ipt], pTbin[ipt+1]) );
                if(isoboth || ipt >= 23) // isoboth or >16GeV don't subtract background
                  FitMinv(h_minv, npion[isoboth][isopair][pttype], enpion[isoboth][isopair][pttype], false, 0.10,0.17);
                else if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
                  FitMinv(h_minv, npion[isoboth][isopair][pttype], enpion[isoboth][isopair][pttype], true, 0.11,0.16);
                else if(ipt < 23)  // <16GeV subtract background
                  FitMinv(h_minv, npion[isoboth][isopair][pttype], enpion[isoboth][isopair][pttype], true, 0.10,0.17);
                if(isoboth)
                  npion[isoboth][isopair][pttype] /= 1.1;
                else                             
                  npion[isoboth][isopair][pttype] /= bck[part/2][ipt] * meff[part/2][ipt];
                delete h_minv;
              }

        double Miss, eMiss, MissEta, eMissEta, Merge1, eMerge1, Merge2, eMerge2, BadPass, eBadPass, Veto, eVeto;
        qt_miss->Query(ipt, part, xpt, Miss, eMiss);
        qt_miss_eta->Query(ipt, part, xpt, MissEta, eMissEta);
        qt_merge1->Query(ipt, part, xpt, Merge1, eMerge1);
        qt_merge2->Query(ipt, part, xpt, Merge2, eMerge2);
        qt_badpass->Query(ipt, part/2, xpt, BadPass, eBadPass);
        qt_veto->Query(ipt, part, xpt, Veto, eVeto);

        double nphoton = nph[0];
        double enphoton = enph[0];
        double n2photon = npion[0][0][0];
        double en2photon = enpion[0][0][0];
        double n2photon2pt = npion[0][0][1];
        double en2photon2pt = enpion[0][0][1];
        double Eff = Conv[part] * Prob * ToF[part];
        double ASee = A * (1+MissEta)/(1+2.*MissEta) * (1+2.*Miss+Merge1);
        ndir[0][part] = nphoton/Eff - (1 + Miss + Merge1*Conv[part]*(1-Conv[part]) + ASee) * n2photon/Eff/Eff - Merge2/2.*BadPass * n2photon2pt/Eff/Eff;
        endir[0][part] = sqrt(pow(eToF[part],2)*pow((2*(1 + (1 - Conv[part])*Conv[part]*Merge1 + Miss + (A*(1 + Merge1 + 2.*Miss)*(1 + MissEta))/(1 + 2.*MissEta))*n2photon)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],3)) + (1.*BadPass*Merge2*n2photon2pt)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],3)) - nphoton/(Conv[part]*Prob*pow(ToF[part],2)),2) + pow(eProb,2)*pow((2*(1 + (1 - Conv[part])*Conv[part]*Merge1 + Miss + (A*(1 + Merge1 + 2.*Miss)*(1 + MissEta))/(1 + 2.*MissEta))*n2photon)/(pow(Conv[part],2)*pow(Prob,3)*pow(ToF[part],2)) + (1.*BadPass*Merge2*n2photon2pt)/(pow(Conv[part],2)*pow(Prob,3)*pow(ToF[part],2)) - nphoton/(Conv[part]*pow(Prob,2)*ToF[part]),2) + pow(eConv[part],2)*pow(-((((1 - Conv[part])*Merge1 - Conv[part]*Merge1)*n2photon)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2))) + (2*(1 + (1 - Conv[part])*Conv[part]*Merge1 + Miss + (A*(1 + Merge1 + 2.*Miss)*(1 + MissEta))/(1 + 2.*MissEta))*n2photon)/(pow(Conv[part],3)*pow(Prob,2)*pow(ToF[part],2)) + (1.*BadPass*Merge2*n2photon2pt)/(pow(Conv[part],3)*pow(Prob,2)*pow(ToF[part],2)) - nphoton/(pow(Conv[part],2)*Prob*ToF[part]),2) + (0.25*pow(BadPass,2)*pow(en2photon2pt,2)*pow(Merge2,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (pow(en2photon,2)*pow(1 + (1 - Conv[part])*Conv[part]*Merge1 + Miss + (A*(1 + Merge1 + 2.*Miss)*(1 + MissEta))/(1 + 2.*MissEta),2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (pow(eA,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(1 + MissEta,2)*pow(n2photon,2))/(pow(Conv[part],4)*pow(1 + 2.*MissEta,2)*pow(Prob,4)*pow(ToF[part],4)) + (pow(eMissEta,2)*pow((-2.*A*(1 + Merge1 + 2.*Miss)*(1 + MissEta))/pow(1 + 2.*MissEta,2) + (A*(1 + Merge1 + 2.*Miss))/(1 + 2.*MissEta),2)*pow(n2photon,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (pow(eMerge1,2)*pow((1 - Conv[part])*Conv[part] + (A*(1 + MissEta))/(1 + 2.*MissEta),2)*pow(n2photon,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (pow(eMiss,2)*pow(1 + (2.*A*(1 + MissEta))/(1 + 2.*MissEta),2)*pow(n2photon,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (0.25*pow(BadPass,2)*pow(eMerge2,2)*pow(n2photon2pt,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (0.25*pow(eBadPass,2)*pow(Merge2,2)*pow(n2photon2pt,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + pow(enphoton,2)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2)));

        nphoton = nph[1];
        enphoton = enph[1];
        double nisoboth = npion[1][0][0];
        double enisoboth = enpion[1][0][0];
        double nisopair = npion[0][1][0];
        double enisopair = enpion[0][1][0];
        double nisopair2pt = npion[0][1][1];
        double enisopair2pt = enpion[0][1][1];
        double AIso = A * Veto * (1.+MissEta)/(1.+2.*MissEta) * (1+2.*Miss+Merge1);
        ndir[1][part] = nphoton/Conv[part] - (1. + Merge1*Conv[part]*(1.-Conv[part])) * nisoboth/pow(Conv[part],2) - Miss * nisopair/pow(Conv[part],2) - Merge2/2.*BadPass * nisopair2pt - AIso * nisopair;
        endir[1][part] = sqrt((pow(enisoboth,2)*pow(1 + (1 - Conv[part])*Conv[part]*Merge1,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (0.25*pow(BadPass,2)*pow(enisopair2pt,2)*pow(Merge2,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (pow(A,2)*pow(eVeto,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(nisopair,2))/(pow(Conv[part],4)*pow(1 + 2.*MissEta,2)*pow(Prob,4)*pow(ToF[part],4)) + (0.25*pow(BadPass,2)*pow(eMerge2,2)*pow(nisopair2pt,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (0.25*pow(eBadPass,2)*pow(Merge2,2)*pow(nisopair2pt,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + pow(enphoton,2)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2)) + (pow(eA,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(nisopair,2)*pow(MissEta + Veto,2))/(pow(Conv[part],4)*pow(1 + 2.*MissEta,2)*pow(Prob,4)*pow(ToF[part],4)) + pow(eToF[part],2)*pow((2*(1 + (1 - Conv[part])*Conv[part]*Merge1)*nisoboth)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],3)) + (2*Miss*nisopair)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],3)) + (1.*BadPass*Merge2*nisopair2pt)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],3)) - nphoton/(Conv[part]*Prob*pow(ToF[part],2)) + (2*A*(1 + Merge1 + 2.*Miss)*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],3)),2) + pow(eProb,2)*pow((2*(1 + (1 - Conv[part])*Conv[part]*Merge1)*nisoboth)/(pow(Conv[part],2)*pow(Prob,3)*pow(ToF[part],2)) + (2*Miss*nisopair)/(pow(Conv[part],2)*pow(Prob,3)*pow(ToF[part],2)) + (1.*BadPass*Merge2*nisopair2pt)/(pow(Conv[part],2)*pow(Prob,3)*pow(ToF[part],2)) - nphoton/(Conv[part]*pow(Prob,2)*ToF[part]) + (2*A*(1 + Merge1 + 2.*Miss)*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,3)*pow(ToF[part],2)),2) + pow(enisopair,2)*pow(-(Miss/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2))) - (A*(1 + Merge1 + 2.*Miss)*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2)),2) + pow(eMissEta,2)*pow(-((A*(1 + Merge1 + 2.*Miss)*nisopair)/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2))) + (2.*A*(1 + Merge1 + 2.*Miss)*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*pow(1 + 2.*MissEta,2)*pow(Prob,2)*pow(ToF[part],2)),2) + pow(eMiss,2)*pow(-(nisopair/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2))) - (2.*A*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2)),2) + pow(eMerge1,2)*pow(-(((1 - Conv[part])*nisoboth)/(Conv[part]*pow(Prob,2)*pow(ToF[part],2))) - (A*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2)),2) + pow(eConv[part],2)*pow(-((((1 - Conv[part])*Merge1 - Conv[part]*Merge1)*nisoboth)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2))) + (2*(1 + (1 - Conv[part])*Conv[part]*Merge1)*nisoboth)/(pow(Conv[part],3)*pow(Prob,2)*pow(ToF[part],2)) + (2*Miss*nisopair)/(pow(Conv[part],3)*pow(Prob,2)*pow(ToF[part],2)) + (1.*BadPass*Merge2*nisopair2pt)/(pow(Conv[part],3)*pow(Prob,2)*pow(ToF[part],2)) - nphoton/(pow(Conv[part],2)*Prob*ToF[part]) + (2*A*(1 + Merge1 + 2.*Miss)*nisopair*(MissEta + Veto))/(pow(Conv[part],3)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2)),2));
      } // part

      for(int isolated=0; isolated<2; isolated++)
      {
        Chi2Fit(3, ndir[isolated], endir[isolated], nsys[isolated][isys], ensys[isolated][isys]);
        if(isys > 0)
        {
          int index = isolated + 2*(isys-1);
          double rsys = nsys[isolated][isys]/nsys[isolated][0];
          double ersys = rsys*sqrt(pow(ensys[isolated][isys]/nsys[isolated][isys],2) + pow(ensys[isolated][0]/nsys[isolated][0],2));
          if( TMath::Finite(rsys+ersys) )
            qt_sys->Fill(ipt, index, xpt, rsys, ersys);
        }
      } // isolated
    } // isys
  } // ipt

  qt_sys->Write();
  for(int isoboth=0; isoboth<2; isoboth++)
    for(int isopair=0; isopair<2; isopair++)
      for(int pttype=0; pttype<2; pttype++)
        if( (1-isoboth)*(1-isopair) == 1 ||
            isoboth*(1-isopair)*(1-pttype) == 1 ||
            (1-isoboth)*isopair == 1 )
          for(int part=0; part<3; part++)
            for(int isys=0; isys<4; isys++)
            {
              int ic = isoboth + 2*isopair + 2*2*pttype + 2*2*2*part + 2*2*2*3*isys;
              mcw( ic, Form("Minv-isoboth%d-isopair%d-pttype%d-part%d-isys%d",
                    isoboth,isopair,pttype,part,isys) );
            }
  qt_sys->Close();
}
