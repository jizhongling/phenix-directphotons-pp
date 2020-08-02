#include "Pileup.h"
#include "QueryTree.h"
#include "DataBase.h"

void anaPileup_IsoPhoton(const int process = 0)
{
  const int nThread = 100;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist-DC3sigma.txt");

  QueryTree *qt_pile = new QueryTree(Form("histos/Pileup-isophoton-%d.root",process), "RECREATE");

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread )
      continue;
    else if( thread >= (process+1)*nThread )
      break;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15811/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("h_events");

    // h[ic][part]
    TH1 *h_1photon[2][3];

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
