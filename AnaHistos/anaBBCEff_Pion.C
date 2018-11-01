#include "GlobalVars.h"
#include "BBCCounts.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void anaBBCEff_Pion(const int process = 0)
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  const int nThread = 20;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  TFile *f_out = new TFile(Form("pileup/BBCEff-pion-%d.root",process), "RECREATE");
  for(int ic=0; ic<2; ic++)
    mc(ic, 6,5);

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

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    THnSparse *hn_trig = (THnSparse*)f->Get("hn_bbc_pion");
    TAxis *axis_sec = hn_trig->GetAxis(0);
    TAxis *axis_pt = hn_trig->GetAxis(1);
    TAxis *axis_minv = hn_trig->GetAxis(2);
    TAxis *axis_cond = hn_trig->GetAxis(3);

    ULong64_t nclock = db->GetClockLive(runnumber);
    ULong64_t nmb = db->GetBBCNovtxLive(runnumber);

    for(int part=0; part<2; part++)
    {
      for(int ipt=0; ipt<npT; ipt++)
      {
        int ig = ipt*2+part;

        axis_sec->SetRange(secl[part],sech[part]);
        axis_pt->SetRange(ipt+1,ipt+1);
        TH1 *h_minv;

        double nt, ent;
        mcd(0, ipt+1);
        axis_cond->SetRange(1,1);
        h_minv = hn_trig->Projection(2);
        h_minv->Rebin(10);
        h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
        FitMinv(h_minv, nt, ent);
        delete h_minv;

        double np, enp;
        mcd(1, ipt+1);
        axis_cond->SetRange(2,2);
        h_minv = hn_trig->Projection(2);
        h_minv->Rebin(10);
        h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
        FitMinv(h_minv, np, enp);
        delete h_minv;

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

      f_out->cd();
      mcw( 0, Form("Run%d-part%d-total",runnumber,part) );
      mcw( 1, Form("Run%d-part%d-passed",runnumber,part) );
      for(int ic=0; ic<2; ic++)
        gROOT->ProcessLine( Form("c%d->Clear(\"D\");",ic) );
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
