#include "Pileup.h"
#include "QueryTree.h"
#include "DataBase.h"
#include "FitMinv.h"

void anaPileup_Isolated(const int process = 0)
{
  const int nThread = 100;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist-DC3sigma.txt");

  QueryTree *qt_pile = new QueryTree(Form("histos/Pileup-isolated-%d.root",process), "RECREATE");

  DataBase *db = new DataBase();

  mc();
  mcd();

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

    // h[evtype][part]
    TH1 *h_1photon[3][3];
    TH2 *h2_isoboth[3][3];
    TH2 *h2_isopair[3][3];

    int evtype = 2;
    int bbc10cm = 1;
    int checkmap = 1;
    int ival = 1;

    TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
    h_1photon_t = (TH1*)h_1photon_t->Clone();
    h_1photon_t->Reset();

    TH2 *h2_2photon_t = (TH2*)f->Get("h2_2photon_0");
    h2_2photon_t = (TH2*)h2_2photon_t->Clone();
    h2_2photon_t->Reset();

    for(int part=0; part<3; part++)
    {
      h_1photon[evtype][part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_type%d_part%d",evtype,part));
      h2_isoboth[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isoboth_type%d_part%d",evtype,part));
      h2_isopair[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isopair_type%d_part%d",evtype,part));

      for(int pattern=0; pattern<3; pattern++)
        for(int evenodd=0; evenodd<2; evenodd++)
        {
          int isolated = 1;
          int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*bbc10cm + 3*2*3*4*2*checkmap + 3*2*3*4*2*2*isolated + 3*2*3*4*2*2*2*ival;
          TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
          h_1photon[evtype][part]->Add(h_tmp);
          delete h_tmp;

          for(int isopair=0; isopair<2; isopair++)
          {
            int isoboth = 1;
            int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*bbc10cm + 3*2*3*4*2*checkmap + 3*2*3*4*2*2*isoboth + 3*2*3*4*2*2*2*isopair + 3*2*3*4*2*2*2*2*ival;
            TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
            h2_isoboth[evtype][part]->Add(h2_tmp);
            delete h2_tmp;
          } // isopair

          for(int isoboth=0; isoboth<2; isoboth++)
          {
            int isopair = 1;
            int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*bbc10cm + 3*2*3*4*2*checkmap + 3*2*3*4*2*2*isoboth + 3*2*3*4*2*2*2*isopair + 3*2*3*4*2*2*2*2*ival;
            TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
            h2_isopair[evtype][part]->Add(h2_tmp);
            delete h2_tmp;
          } // isoboth
        }
    }

    unsigned long long nclock = db->GetClockLive(runnumber);
    unsigned long long nmb = db->GetBBCNarrowLive(runnumber);
    unsigned long long scaledown = db->GetERT4x4cScaledown(runnumber) + 1;

    double nev = h_events->GetBinContent( h_events->GetXaxis()->FindBin("ert_c_10cm") );
    //double nev = nmb / scaledown;

    for(int ipt=0; ipt<npT; ipt++)
      for(int part=0; part<3; part++)
      {
        double xx = (double)nmb / (double)nclock;
        double yy[3], eyy[3];

        yy[0] = h_1photon[evtype][part]->IntegralAndError(pTlow[0][ipt],pThigh[0][ipt], eyy[0]);

        TH1 *h_minv;

        h_minv = (TH1*)h2_isoboth[evtype][part]->ProjectionY("h_py", pTlow[0][ipt],pThigh[0][ipt])->Clone("h_minv");
        h_minv->Rebin(10);
        FitMinv(h_minv, yy[1], eyy[1], false, 0.10,0.17);
        delete h_minv;

        h_minv = (TH1*)h2_isopair[evtype][part]->ProjectionY("h_py", pTlow[0][ipt],pThigh[0][ipt])->Clone("h_minv");
        h_minv->Rebin(10);
        FitMinv(h_minv, yy[2], eyy[2], false, 0.10,0.17);
        delete h_minv;

        for(int iph=0; iph<3; iph++)
          if( TMath::Finite(yy[iph]+eyy[iph]) )
          {
            int ig = iph + 3*part + 3*3*ipt;
            qt_pile->Fill(runnumber, ig, xx, yy[iph]/nev, eyy[iph]/nev);
          }
      }

    delete f;
  }

  qt_pile->Save();
}
