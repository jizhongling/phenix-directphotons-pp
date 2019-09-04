#include "DataBase.h"
#include "FitMinv.h"

void anaPileup_Sasha(const int process = 0)
{
  const int nThread = 20;
  int thread = -1;
  int irun = 0;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  TGraphErrors *gr[4];
  int igp[4] = {};
  for(int ig=0; ig<4; ig++)
  {
    mc(ig, 5,4);
    gr[ig] = new TGraphErrors(nThread);
    gr[ig]->SetName(Form("gr_%d",ig));
  }

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510MinBias/12233/data/Pi0PP-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("hevtype");
    TH1 *h_minv[2][2];  // h_minv[ic][is]
    double npion[2][2];  // npion[ic][is]
    double enpion[2][2];  // enpion[ic][is]
    for(int ic=0; ic<2; ic++)
      for(int is=0; is<2; is++)
      {
        h_minv[ic][is] = new TH1F(Form("h_minv_%d",ic*2+is), Form("#%d",runnumber), 1000,0.,1.);
        npion[ic][is] = 0.;
        enpion[ic][is] = 0.;
      }

    const char *cname[2] = {"p", "tp"};
    for(int ic=0; ic<2; ic++)
      for(int is=0; is<3; is++)
        for(int ip=4; ip<20; ip++)
        {
          TH1 *mchist = (TH1*)f->Get(Form("mc_s%d_bcc0_pt_%03d_%s",is,5*ip,cname[ic]));
          h_minv[ic][is/2]->Add(mchist);
        }

    TF1 *fn_fit = new TF1("fn_fit", "gaus+pol2(3)", 0.06, 0.25);
    TF1 *fn_bg = new TF1("fn_bg", "pol2", 0.06, 0.25);

    unsigned long long nclock = db->GetClockLive(runnumber);
    unsigned long long nmb = db->GetBBCNarrowLive(runnumber);
    double nev = h_events->GetBinContent( h_events->GetXaxis()->FindBin("bbc_novtx_10cm") );

    for(int ic=0; ic<2; ic++)
      for(int is=0; is<2; is++)
      {
        int ig = ic*2+is;
        mcd(ig, irun+1);
        h_minv[ic][is]->Rebin(10);
        h_minv[ic][is]->Scale(0.5);
        FitMinv(h_minv[ic][is], npion[ic][is], enpion[ic][is]);
      }

    for(int ic=0; ic<2; ic++)
      for(int is=0; is<2; is++)
      {
        double xx = (double)nmb/(double)nclock;
        double yy = npion[ic][is] / nev;
        double eyy = enpion[ic][is] / nev;
        if( TMath::Finite(yy+eyy) )
        {
          gr[ig]->SetPoint(igp[ig], xx, yy);
          gr[ig]->SetPointError(igp[ig], 0., eyy);
          igp[ig]++;
        }
        delete h_minv[ic][is];
      }
    delete f;
    irun++;
  }

  TFile *f_out = new TFile(Form("histos/Sasha-%d.root",process), "RECREATE");
  for(int ig=0; ig<4; ig++)
  {
    gr[ig]->Set(igp[ig]);
    mcw( ig, Form("proc%d-cond%d", process, ig) );
    gr[ig]->Write();
  }
  f_out->Close();
}
