#include "GlobalVars.h"
#include "BBCCounts.h"
#include "GetEfficiency.h"

void anaBBCEff_Photon(const Int_t process = 0)
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  const Int_t nThread = 1000;
  Int_t thread = -1;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  TFile *f_out = new TFile(Form("pileup/BBCEff-photon-%d.root",process), "RECREATE");

  TGraphAsymmErrors *gr[npT*2];
  Int_t igp[npT*2] = {};
  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      Int_t ig = ipt*2+part;
      gr[ig] = new TGraphAsymmErrors(nThread);
      gr[ig]->SetName(Form("gr_%d",ig));
    }

  ReadClockCounts();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;
    if( thread%10 == 0 ) cout << "nfile = " << thread << endl;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH3 *h3_trig = (TH3*)f->Get("h3_bbc");

    ULong64_t nclock = GetClockLive(runnumber);
    ULong64_t nmb = GetBBCNovtxLive(runnumber);

    for(Int_t part=0; part<2; part++)
      for(Int_t ipt=0; ipt<npT; ipt++)
      {
        Int_t ig = ipt*2+part;

        Double_t nt = h3_trig->Integral(secl[part],sech[part], ipt+1,ipt+1, 1,1);
        Double_t np = h3_trig->Integral(secl[part],sech[part], ipt+1,ipt+1, 2,2);

        Double_t xx = (Double_t)nmb / (Double_t)nclock;
        Double_t yy, eyyl, eyyh;
        if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
        {
          eyyl = yy * sqrt( pow(ent/nt,2.) + pow(enp/np,2.) );
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
  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      Int_t ig = ipt*2+part;
      gr[ig]->Set(igp[ig]);
      gr[ig]->Write();
    }
  f_out->Close();
}
