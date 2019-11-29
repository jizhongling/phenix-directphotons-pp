#include "Pileup.h"
#include "QueryTree.h"
#include "DataBase.h"
#include "FitMinv.h"

void anaPileup(const int process = 0)
{
  const int nThread = 20;
  int thread = -1;
  int irun = 0;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  QueryTree *qt_pile = new QueryTree(Form("histos/Pileup-%d.root",process), "RECREATE");

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread )
      continue;
    else if( thread >= (process+1)*nThread )
      break;

    TFile *f_ert = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    TFile *f_mb = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-MB/PhotonNode-%d.root",runnumber));
    if( f_ert->IsZombie() || f_mb->IsZombie() ) continue;

    TH1 *h_events_ert = (TH1*)f_ert->Get("h_events");
    TH1 *h_events_mb = (TH1*)f_mb->Get("h_events");
    THnSparse *hn_pion[2];
    hn_pion[0] = (THnSparse*)f_ert->Get("hn_pion");
    hn_pion[1] = (THnSparse*)f_mb->Get("hn_pion");
    hn_pion[0]->GetAxis(4)->SetRange(3,3);
    hn_pion[1]->GetAxis(4)->SetRange(1,1);

    unsigned long long nclock = db->GetClockLive(runnumber);
    unsigned long long nmb = db->GetBBCNarrowLive(runnumber);
    unsigned long long scaledown = db->GetERT4x4cScaledown(runnumber) + 1;
    
    double nev[2];
    //nev[0] = h_events_ert->GetBinContent( h_events_ert->GetXaxis()->FindBin("ert_c") );
    nev[0] = nmb / scaledown;
    nev[1] = h_events_mb->GetBinContent( h_events_mb->GetXaxis()->FindBin("bbc_narrow_10cm") );

    TF1 *fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", 0.06, 0.25);
    TF1 *fn_bg = new TF1("fn_bg", "pol2", 0.06, 0.25);

    for(int ipt=0; ipt<npT; ipt++)
      for(int id=0; id<2; id++)
        for(int ic=0; ic<2; ic++)
          for(int is=0; is<2; is++)
          {
            int ig = ipt*8+id*4+ic*2+is;
            mcd(ig, irun+1);
            
            double npion, enpion;

            hn_pion[id]->GetAxis(3)->SetRange(ic+3,ic+3);
            hn_pion[id]->GetAxis(0)->SetRange(secl[is],sech[is]);
            hn_pion[id]->GetAxis(1)->SetRange(pTlow[id][ipt], pThigh[id][ipt]);

            TH1 *h_minv = hn_pion[id]->Projection(2);
            h_minv->Rebin(10);
            h_minv->SetTitle(Form("#%d",runnumber));
            FitMinv(h_minv, npion, enpion);
            delete h_minv;

            double xx = (double)nmb / (double)nclock;
            double yy = npion / nev[id];
            double eyy = enpion / nev[id];
            if( TMath::Finite(yy+eyy) )
              qt_pile->Fill(runnumber, ig, xx, yy, eyy);
          }

    delete f_ert;
    delete f_mb;
    irun++;
  }

  qt_pile->Write();
  for(int ipt=0; ipt<npT; ipt++)
    for(int id=0; id<2; id++)
      for(int ic=0; ic<2; ic++)
        for(int is=0; is<2; is++)
          mcw( ig, Form("proc%d-data%d-cond%d-pt%d-%d", process, id, ic*2+is, pTlow[id][ipt], pThigh[id][ipt]) );
  qt_pile->Close();
}
