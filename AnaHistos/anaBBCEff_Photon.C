#include "GlobalVars.h"
#include "DataBase.h"
#include "GetEfficiency.h"

void anaBBCEff_Photon(const int process = 0)
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  const int nThread = 1000;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  TFile *f_out = new TFile(Form("pileup/BBCEff-photon-%d.root",process), "RECREATE");

  TGraphAsymmErrors *gr[npT*2];
  int igp[npT*2] = {};
  for(int part=0; part<2; part++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      int ig = ipt*2+part;
      gr[ig] = new TGraphAsymmErrors(nThread);
      gr[ig]->SetName(Form("gr_%d",ig));
    }

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;
    if( thread%10 == 0 ) cout << "nfile = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH3 *h3_trig = (TH3*)f->Get("h3_bbc");

    ULong64_t nclock = db->GetClockLive(runnumber);
    ULong64_t nmb = db->GetBBCNovtxLive(runnumber);

    for(int part=0; part<2; part++)
      for(int ipt=0; ipt<npT; ipt++)
      {
        int ig = ipt*2+part;

        double nt = h3_trig->Integral(secl[part],sech[part], ipt+1,ipt+1, 1,1);
        double np = h3_trig->Integral(secl[part],sech[part], ipt+1,ipt+1, 2,2);

        double xx = (double)nmb / (double)nclock;
        double yy, eyyl, eyyh;
        if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
        {
          eyyl = yy * sqrt( pow(ent/nt,2) + pow(enp/np,2) );
          eyyh = 0.;
        }
        if( yy >= 0. && eyyl >= 0. && eyyl < TMath::Infinity() )
        {
          gr[ig]->SetPoint(igp[ig], xx, yy);
          gr[ig]->SetPointError(igp[ig], 0.,0., eyyl,eyyh);
          igp[ig]++;
        }
      }

    delete f;
  }

  f_out->cd();
  for(int part=0; part<2; part++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      int ig = ipt*2+part;
      gr[ig]->Set(igp[ig]);
      gr[ig]->Write();
    }
  f_out->Close();
}
