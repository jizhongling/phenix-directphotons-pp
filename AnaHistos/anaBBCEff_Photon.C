#include "GlobalVars.h"
#include "QueryTree.h"
#include "DataBase.h"
#include "GetEfficiency.h"

void anaBBCEff_Photon(const int process = 0)
{
  const int nThread = 1000;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  QueryTree *qt_bbc = new QueryTree(Form("histos/BBCEff-photon-%d.root",process), "RECREATE");

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread )
      continue;
    else if( thread >= (process+1)*nThread )
      break;
    if( thread && thread%10 == 0 )
      cout << "nfile = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH3 *h3_trig = (TH3*)f->Get("h3_bbc");

    unsigned long long nclock = db->GetClockLive(runnumber);
    unsigned long long nmb = db->GetBBCNovtxLive(runnumber);

    for(int part=0; part<2; part++)
      for(int ipt=0; ipt<npT; ipt++)
      {
        int ig = ipt*2+part;

        double nt = h3_trig->Integral(secl[part],sech[part], ipt+1,ipt+1, 1,1);
        double np = h3_trig->Integral(secl[part],sech[part], ipt+1,ipt+1, 2,2);

        double xpt = (double)nmb / (double)nclock;
        double yy, eyyl, eyyh;
        if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
        {
          eyyl = yy * sqrt( pow(ent/nt,2) + pow(enp/np,2) );
          eyyh = 0.;
        }
        if( TMath::Finite(yy+eyyl+eyyh) )
          qt_bbc->Fill(runnumber, ig, xpt, yy, eyyl, eyyh);
      }

    delete f;
  }

  qt_bbc->Save();
}
