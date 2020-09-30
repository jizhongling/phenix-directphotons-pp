#include "Pileup.h"
#include "QueryTree.h"
#include "DataBase.h"
#include "FitMinv.h"

void anaPileup(const int process = 0)
{
  const int nThread = 100;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist-DC3sigma.txt");

  QueryTree *qt_pile = new QueryTree(Form("histos/Pileup-%d.root",process), "RECREATE");

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

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/16669/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("h_events");

    // h[iso][part]
    TH1 *h_1photon[2][3];
    TH2 *h2_2photon[3][3];

    int evtype = 2;
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
      for(int iso=0; iso<3; iso++)
      {
        if(iso<2)
          h_1photon[iso][part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_iso%d_part%d",iso,part));
        h2_2photon[iso][part] = (TH2*)h2_2photon_t->Clone(Form("h2_2photon_iso%d_part%d",iso,part));
      }

      for(int isolated=0; isolated<2; isolated++)
      {
        int ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isolated + 3*3*2*2*ival;
        TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
        h_1photon[0][part]->Add(h_tmp);
        if(isolated==1)
          h_1photon[1][part]->Add(h_tmp);
        delete h_tmp;
      }

      for(int isoboth=0; isoboth<2; isoboth++)
        for(int isopair=0; isopair<2; isopair++)
        {
          int ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isoboth + 3*3*2*2*isopair + 3*3*2*2*2*ival;
          TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
          h2_2photon[0][part]->Add(h2_tmp);
          if(isoboth==1)
            h2_2photon[1][part]->Add(h2_tmp);
          if(isopair==1)
            h2_2photon[2][part]->Add(h2_tmp);
          delete h2_tmp;
        } // isoboth, isopair
    }

    unsigned long long nclock = db->GetClockLive(runnumber);
    unsigned long long nmb = db->GetBBCNarrowLive(runnumber);
    unsigned long long scaledown = db->GetERT4x4cScaledown(runnumber) + 1;

    double nev = h_events->GetBinContent( h_events->GetXaxis()->FindBin("ert_c_10cm") );
    //double nev = nmb / scaledown;

    for(int part=0; part<3; part++)
      for(int ipt=0; ipt<npT; ipt++)
      {
        double xx = (double)nmb / (double)nclock;
        double yy[5] = {}, eyy[5] = {};

        for(int iso=0; iso<3; iso++)
        {
          TH1 *h_minv;
          h_minv = (TH1*)h2_2photon[iso][part]->ProjectionY("h_py", pTlow[0][ipt],pThigh[0][ipt])->Clone("h_minv");
          h_minv->Rebin(10);
          FitMinv(h_minv, yy[2+iso], eyy[2+iso], false, 0.10,0.17);
          delete h_minv;

          if(iso<2)
          {
            yy[iso] = h_1photon[iso][part]->IntegralAndError(pTlow[0][ipt],pThigh[0][ipt], eyy[iso]);
            yy[iso] -= yy[2+iso];
            eyy[iso] = sqrt(eyy[iso]*eyy[iso] + eyy[iso+2]*eyy[iso+2]);
          }
        } // iso

        for(int iph=0; iph<5; iph++)
          if( TMath::Finite(yy[iph]+eyy[iph]) && eyy[iph] > 0. )
          {
            int ig = iph + 5*part + 5*3*ipt;
            qt_pile->Fill(runnumber, ig, xx, yy[iph]/nev, eyy[iph]/nev);
          } // iph
      } // part, ipt

    delete f;
  }

  qt_pile->Save();
}
