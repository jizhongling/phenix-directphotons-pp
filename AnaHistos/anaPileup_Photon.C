#include "Pileup.h"
#include "QueryTree.h"
#include "DataBase.h"

void anaPileup_Photon(const int process = 0)
{
  const int nThread = 100;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  const double A = 0.28;
  const double eA = 0.05;

  QueryTree *qt_pile = new QueryTree(Form("histos/Pileup-photon-%d.root",process), "RECREATE");

  QueryTree *qt_misscorr = new QueryTree("data/MissCorr.root");
  QueryTree *qt_merge1 = new QueryTree("data/MergeCorr-1photon.root");
  QueryTree *qt_merge2 = new QueryTree("data/MergeCorr-2photon.root");

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15410/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("h_events");

    // h[part]
    TH1 *h_1photon[3];
    TH2 *h2_2photon[3];

    int bbc10cm = 1;
    int evtype = 2;
    int ival = 1;

    TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
    h_1photon_t->Reset();
    for(int part=0; part<3; part++)
    {
      h_1photon[part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_part%d",part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int evenodd=0; evenodd<2; evenodd++)
          for(int pattern=0; pattern<3; pattern++)
            for(int isolated=0; isolated<2; isolated++)
            {
              int ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isolated + 8*2*3*2*evtype + 8*2*3*2*3*bbc10cm + 8*2*3*2*3*2*ival + 8*2*3*2*3*2*4*(tof-1);
              TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
              h_1photon[part]->Add(h_tmp);
              delete h_tmp;
            }
    }

    TH2 *h2_2photon_t = (TH2*)f->Get("h2_2photon_0");
    h2_2photon_t->Reset();
    for(int part=0; part<3; part++)
    {
      h2_2photon[part] = (TH2*)h2_2photon_t->Clone(Form("h2_2photon_part%d",part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int evenodd=0; evenodd<2; evenodd++)
          for(int pattern=0; pattern<3; pattern++)
            for(int isoboth=0; isoboth<2; isoboth++)
              for(int isopair=0; isopair<2; isopair++)
              {
                int ih = sector + 8*evenodd + 8*2*pattern + 8*2*3*isoboth + 8*2*3*2*isopair + 8*2*3*2*2*evtype + 8*2*3*2*2*3*bbc10cm + 8*2*3*2*2*3*2*ival + 8*2*3*2*2*3*2*4*(tof-1);
                TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
                h2_2photon[part]->Add(h2_tmp);
                delete h2_tmp;
              }
    }

    unsigned long long nclock = db->GetClockLive(runnumber);
    unsigned long long nmb = db->GetBBCNarrowLive(runnumber);
    unsigned long long scaledown = db->GetERT4x4cScaledown(runnumber) + 1;

    //double nev = h_events->GetBinContent( h_events->GetXaxis()->FindBin("ert_c_10cm") );
    double nev = nmb / scaledown;

    int id = 0;
    int ic = 1;
    for(int ipt=0; ipt<npT; ipt++)
      for(int part=0; part<3; part++)
      {
        int ig = part + 3*ic + 2*3*id + 2*2*3*ipt;

        double nphoton = h_1photon[part]->Integral(pTlow[id][ipt], pThigh[id][ipt]);
        double npion = h2_2photon[part]->ProjectionY("h_minv", pTlow[id][ipt], pThigh[id][ipt])->Integral(111,160);
        npion /= 1.1;

        double xpt, MissCorr, eMissCorr, MergeCorr1, eMergeCorr1, MergeCorr2, eMergeCorr2, Veto, eVeto;
        qt_misscorr->Query(ipt, part, xpt, MissCorr, eMissCorr);
        qt_mergecorr1->Query(ipt, part, xpt, MergeCorr1, eMergeCorr1);
        qt_mergecorr2->Query(ipt, part, xpt, MergeCorr2, eMergeCorr2);

        double ndir = nphoton - MissPass * npion - A * ( 1. + MissAll ) * npion;
        double endir = sqrt( nphoton + pow(MissPass + A * ( 1. + MissAll ), 2) * npion );

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
