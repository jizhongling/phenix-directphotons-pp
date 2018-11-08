#include "GlobalVars.h"
#include "DataBase.h"
#include "GetEfficiency.h"

void anaToFEff_Photon(const int process = 0)
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  const int nThread = 1000;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  TFile *f_out = new TFile(Form("pileup/ToF-photon-%d.root",process), "RECREATE");

  TGraphAsymmErrors *gr[2];
  int igp[2] = {};
  for(int part=0; part<2; part++)
  {
    gr[part] = new TGraphAsymmErrors(nThread);
    gr[part]->SetName(Form("gr_%d",part));
  }

  DataBase *db = new DataBase();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;
    if( thread%10 == 0 ) cout << "nfile = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
    TAxis *axis_sec = hn_1photon->GetAxis(0);
    TAxis *axis_pt = hn_1photon->GetAxis(1);
    TAxis *axis_cut = hn_1photon->GetAxis(3);
    TAxis *axis_type = hn_1photon->GetAxis(4);

    axis_type->SetRange(3,3);
    axis_pt->SetRange(26,30);

    ULong64_t nclock = db->GetClockLive(runnumber);
    ULong64_t nmb = db->GetBBCNovtxLive(runnumber);

    for(int part=0; part<2; part++)
    {
      axis_sec->SetRange(secl[part],sech[part]);
      TH1 *h_tof = hn_1photon->Projection(3);
      double nt = h_tof->GetBinContent(3);
      double np = h_tof->GetBinContent(4);
      delete h_tof;

      //double xx = (double)nmb / (double)nclock;
      double xx = nt;
      double yy, eyyl, eyyh;
      if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
      {
        eyyl = yy * sqrt( pow(ent/nt,2) + pow(enp/np,2) );
        eyyh = 0.;
      }
      if( !TMath::IsNaN(yy) && !TMath::IsNaN(eyyl) && eyyl >= 0. )
      {
        gr[part]->SetPoint(igp[part], xx, yy);
        gr[part]->SetPointError(igp[part], 0.,0., eyyl,eyyh);
        igp[part]++;
      }
    }

    delete f;
  }

  f_out->cd();
  for(int part=0; part<2; part++)
  {
    gr[part]->Set(igp[part]);
    gr[part]->Write();
  }
  f_out->Close();
}
