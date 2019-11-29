#include "Pileup.h"
#include "QueryTree.h"
#include "DataBase.h"

void anaPileup_IsoPhoton(const int process = 0)
{
  const int nThread = 100;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist-DC3sigma.txt");

  //const double Conv[3] = {0.849, 0.959, 0.959};
  //const double eConv[3] = {0.027, 0.023, 0.023};
  //const double A = 0.28;
  //const double eA = 0.05;

  QueryTree *qt_pile = new QueryTree(Form("histos/Pileup-isophoton-%d.root",process), "RECREATE");

  //QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");
  //QueryTree *qt_miss_eta = new QueryTree("data/MissingRatio-eta.root");
  //QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root");
  //QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root");
  //QueryTree *qt_badpass = new QueryTree("data/MergePassRate.root");
  //QueryTree *qt_veto = new QueryTree("data/SelfVeto.root");

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread )
      continue;
    else if( thread >= (process+1)*nThread )
      break;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15410/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("h_events");

    // h[ic][part]
    TH1 *h_1photon[2][3];
    //TH2 *h2_isoboth[2][3];
    //TH2 *h2_isopair[2][3];
    //TH2 *h2_isopair2pt[2][3];

    int bbc10cm = 1;
    int evtype = 2;
    int ic = 1;
    int ival = 1;

    TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
    h_1photon_t->Reset();
    for(int part=0; part<3; part++)
    {
      h_1photon[ic][part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_ic%d_part%d",ic,part));
      for(int evenodd=0; evenodd<2; evenodd++)
        for(int pattern=0; pattern<3; pattern++)
        {
          int isolated = 1;
          int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*bbc10cm + 3*2*3*4*2*isolated + 3*2*3*4*2*2*ival;
          TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
          h_1photon[ic][part]->Add(h_tmp);
          delete h_tmp;
        }
    }

    //TH2 *h2_2photon_t = (TH2*)f->Get("h2_2photon_0");
    //h2_2photon_t->Reset();
    //for(int part=0; part<3; part++)
    //{
    //  h2_isoboth[ic][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isoboth_ic%d_part%d",ic,part));
    //  for(int evenodd=0; evenodd<2; evenodd++)
    //    for(int pattern=0; pattern<3; pattern++)
    //      for(int isopair=0; isopair<2; isopair++)
    //      {
    //        int isoboth = 1;
    //        int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*bbc10cm + 3*2*3*4*2*isoboth + 3*2*3*4*2*2*isopair + 3*2*3*4*2*2*2*ival;
    //        TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
    //        h2_isoboth[ic][part]->Add(h2_tmp);
    //        delete h2_tmp;
    //      }
    //}

    //for(int part=0; part<3; part++)
    //{
    //  h2_isopair[ic][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isopair_ic%d_part%d",ic,part));
    //  for(int pattern=0; pattern<3; pattern++)
    //    for(int isoboth=0; isoboth<2; isoboth++)
    //    {
    //      int isopair = 1;
    //      int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*bbc10cm + 3*2*3*4*2*isoboth + 3*2*3*4*2*2*isopair + 3*2*3*4*2*2*2*ival;
    //      TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
    //      h2_isopair[ic][part]->Add(h2_tmp);
    //      delete h2_tmp;
    //    }
    //}

    //for(int part=0; part<3; part++)
    //{
    //  h2_isopair2pt[ic][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isopair2pt_ic%d_part%d",ic,part));
    //  for(int pattern=0; pattern<3; pattern++)
    //    for(int isoboth=0; isoboth<2; isoboth++)
    //    {
    //      int isopair = 1;
    //      int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*bbc10cm + 3*2*3*4*2*isoboth + 3*2*3*4*2*2*isopair + 3*2*3*4*2*2*2*ival;
    //      TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
    //      h2_isopair2pt[ic][part]->Add(h2_tmp);
    //      delete h2_tmp;
    //    }
    //}

    unsigned long long nclock = db->GetClockLive(runnumber);
    unsigned long long nmb = db->GetBBCNarrowLive(runnumber);
    unsigned long long scaledown = db->GetERT4x4cScaledown(runnumber) + 1;

    double nev = h_events->GetBinContent( h_events->GetXaxis()->FindBin("ert_c_10cm") );
    //double nev = nmb / scaledown;

    int id = 0;
    for(int ipt=0; ipt<npT; ipt++)
      for(int part=0; part<3; part++)
      {
        int ig = part + 3*ic + 2*3*id + 2*2*3*ipt;

        double nphoton = h_1photon[ic][part]->Integral(pTlow[id][ipt], pThigh[id][ipt]);
        double enphoton = sqrt(nphoton);
        //double nisoboth = h2_isoboth[ic][part]->ProjectionY("h_minv_isoboth", pTlow[id][ipt], pThigh[id][ipt])->Integral(111,160);
        //double nisopair = h2_isopair[ic][part]->ProjectionY("h_minv_isopair", pTlow[id][ipt], pThigh[id][ipt])->Integral(111,160);
        //double nisopair2pt = h2_isopair2pt[ic][part]->ProjectionY("h_minv_isopair2pt", pTlow[id][ipt], pThigh[id][ipt])->Integral(111,160);
        //double enisoboth = sqrt(nisoboth);
        //double enisopair = sqrt(nisopair);
        //double enisopair2pt = sqrt(nisopair2pt);
        //nisoboth /= 1.1;
        //nisopair /= 1.1;
        //nisopair2pt /= 1.1;

        //double xpt, Miss, eMiss, MissEta, eMissEta, Merge1, eMerge1, Merge2, eMerge2, BadPass, eBadPass, Veto, eVeto;
        //qt_miss->Query(ipt, part, xpt, Miss, eMiss);
        //qt_miss_eta->Query(ipt, part, xpt, MissEta, eMissEta);
        //qt_merge1->Query(ipt, part, xpt, Merge1, eMerge1);
        //qt_merge2->Query(ipt, part, xpt, Merge2, eMerge2);
        //qt_badpass->Query(ipt, part/2, xpt, BadPass, eBadPass);
        //qt_veto->Query(ipt, part, xpt, Veto, eVeto);

        //double AIso = A * Veto * (1.+MissEta)/(1.+2.*MissEta) * (1+2.*Miss+Merge1);
        //double ndir = nphoton/Conv[part] - (1. + Merge1*Conv[part]*(1.-Conv[part])) * nisoboth/pow(Conv[part],2) - Miss * nisopair/pow(Conv[part],2) - Merge2/2.*BadPass * nisopair2pt - AIso * nisopair;
        //double endir = sqrt(pow(enphoton,2)/pow(Conv[part],2) + (pow(enisoboth,2)* pow(1. + (1. - Conv[part])*Conv[part]*Merge1,2))/ pow(Conv[part],4) + 0.25*pow(BadPass,2)*pow(enisopair2pt,2)* pow(Merge2,2) + 0.25*pow(BadPass,2)*pow(eMerge2,2)* pow(nisopair2pt,2) + 0.25*pow(eBadPass,2)*pow(Merge2,2)* pow(nisopair2pt,2) + (pow(A,2)*pow(eVeto,2)* pow(1. + Miss,2)* pow(1 + Merge1 + 2.*Miss,2)* pow(nisopair,2))/pow(1. + 2.*Miss,2) + pow(eConv[part],2)* pow(-((((1. - Conv[part])*Merge1 - Conv[part]*Merge1)* nisoboth)/pow(Conv[part],2)) + (2*(1. + (1. - Conv[part])*Conv[part]*Merge1)* nisoboth)/pow(Conv[part],3) + (2*Miss*nisopair)/pow(Conv[part],3) - nphoton/pow(Conv[part],2),2) + (pow(eA,2)*pow(1. + Miss,2)* pow(1 + Merge1 + 2.*Miss,2)* pow(nisopair,2)*pow(Veto,2))/ pow(1. + 2.*Miss,2) + pow(enisopair,2)* pow(-(Miss/pow(Conv[part],2)) - (A*(1. + Miss)*(1 + Merge1 + 2.*Miss)* Veto)/(1. + 2.*Miss),2) + pow(eMerge1,2)* pow(-(((1. - Conv[part])*nisoboth)/Conv[part]) - (A*(1. + Miss)*nisopair*Veto)/ (1. + 2.*Miss),2) + pow(eMiss,2)*pow(-(nisopair/ pow(Conv[part],2)) - (2.*A*(1. + Miss)*nisopair*Veto)/ (1. + 2.*Miss) + (2.*A*(1. + Miss)*(1 + Merge1 + 2.*Miss)* nisopair*Veto)/pow(1. + 2.*Miss,2) - (A*(1 + Merge1 + 2.*Miss)*nisopair*Veto)/ (1. + 2.*Miss),2));

        double ndir = nphoton;
        double endir = enphoton;
        double xx = (double)nmb / (double)nclock;
        double yy = ndir / nev;
        double eyy = endir / nev;
        if( TMath::Finite(yy+eyy) )
          qt_pile->Fill(runnumber, ig, xx, yy, eyy);
      } // ipt, id, ic, part

    delete f;
  }

  qt_pile->Save();
}
