#include "Pileup.h"
#include "BBCCounts.h"
#include "ReadGraph.h"
#include "Chi2Fit.h"
#include "FitMinv.h"

/* Get ipt for TGraph Xaxis gx */
Int_t Get_ipt(Double_t *gx, Double_t xx)
{
  for(Int_t ipt=0; ipt<30; ipt++)
    if( TMath::Abs(gx[ipt] - xx) < 0.2 )
      return ipt;

  cout << "Warning: No matching for pT = " << xx << ", 0 returned!" << endl;
  return 0;
}

void anaPileup_Photon(const Int_t process = 0)
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  const Int_t nThread = 20;
  Int_t thread = -1;
  Int_t irun =0;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  TGraphErrors *gr[npT*8];
  TGraphErrors *gr_run[npT*8];
  Int_t igp[npT*8] = {};
  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t id=0; id<2; id++)
      for(Int_t ic=0; ic<2; ic++)
        for(Int_t is=0; is<2; is++)
        {
          Int_t ig = ipt*8+id*4+ic*2+is;
          mc(ig, 5,4);
          gr[ig] = new TGraphErrors(nThread);
          gr_run[ig] = new TGraphErrors(nThread);
          gr[ig]->SetName(Form("gr_%d",ig));
          gr_run[ig]->SetName(Form("gr_run_%d",ig));
        }

  Double_t xHadron[3][30] = {}, Hadron[3][30] = {}, eHadron[3][30] = {};
  for(Int_t part=0; part<3; part++)
    ReadGraph<TGraphErrors>("data/HadronRatio.root", part, xHadron[part], Hadron[part], eHadron[part]);

  ReadClockCounts();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f_ert = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    TFile *f_mb = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-MB/PhotonNode-%d.root",runnumber));
    if( f_ert->IsZombie() || f_mb->IsZombie() ) continue;

    TH1 *h_events_ert = (TH1*)f_ert->Get("h_events");
    TH1 *h_events_mb = (TH1*)f_mb->Get("h_events");
    THnSparse *hn_1photon[2];
    THnSparse *hn_2photon[2];
    hn_1photon[0] = (THnSparse*)f_ert->Get("hn_1photon");
    hn_1photon[1] = (THnSparse*)f_mb->Get("hn_1photon");
    hn_2photon[0] = (THnSparse*)f_ert->Get("hn_2photon");
    hn_2photon[1] = (THnSparse*)f_mb->Get("hn_2photon");
    hn_1photon[0]->GetAxis(4)->SetRange(3,3);
    hn_1photon[1]->GetAxis(4)->SetRange(1,1);
    hn_2photon[0]->GetAxis(5)->SetRange(3,3);
    hn_2photon[1]->GetAxis(5)->SetRange(1,1);

    ULong64_t nclock = GetClockLive(runnumber);
    ULong64_t nmb = GetBBCNarrowLive(runnumber);
    ULong_t scaledown = GetERT4x4cScaledown(runnumber) + 1;

    Double_t nev[2];
    //nev[0] = h_events_ert->GetBinContent( h_events_ert->GetXaxis()->FindBin("ert_c") );
    nev[0] = nmb / scaledown;
    nev[1] = h_events_mb->GetBinContent( h_events_mb->GetXaxis()->FindBin("bbc_narrow_10cm") );

    TF1 *fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", 0.06, 0.25);
    TF1 *fn_bg = new TF1("fn_bg", "pol2", 0.06, 0.25);

    for(Int_t ipt=0; ipt<npT; ipt++)
      for(Int_t id=0; id<2; id++)
      {
        Double_t Hbar[3] = {}, eHbar[3] = {};
        for(Int_t part=0; part<3; part++)
        {
          Int_t iptL = Get_ipt(xHadron[part], pTbinL[id][ipt]+0.25);
          Chi2Fit(pThigh[id][ipt]-pTlow[id][ipt], &Hadron[part][iptL], &eHadron[part][iptL], Hbar[part], eHbar[part]);
        }
        for(Int_t ic=0; ic<2; ic++)
          for(Int_t is=0; is<2; is++)
          {
            Int_t ig = ipt*8+id*4+ic*2+is;
            mcd(ig, irun+1);

            Double_t nphoton, npion, enpion;

            hn_1photon[id]->GetAxis(3)->SetRange(ic+3,ic+3);
            hn_1photon[id]->GetAxis(1)->SetRange(pTlow[id][ipt], pThigh[id][ipt]);

            TH1 *h_photon = hn_1photon[id]->Projection(0);
            nphoton = h_photon->Integral(secl[is],sech[is]);
            delete h_photon;

            hn_2photon[id]->GetAxis(4)->SetRange(ic+3,ic+3);
            hn_2photon[id]->GetAxis(0)->SetRange(secl[is],sech[is]);
            hn_2photon[id]->GetAxis(1)->SetRange(pTlow[id][ipt], pThigh[id][ipt]);

            TH1 *h_minv = hn_2photon[id]->Projection(2);
            h_minv->Rebin(10);
            h_minv->Scale(0.5);
            h_minv->SetTitle(Form("#%d",runnumber));
            FitMinv(h_minv, npion, enpion);
            delete h_minv;

            Double_t Miss, eMiss;
            if(is==0)
            {
              Chi2Fit(2, Hbar, eHbar, Miss, eMiss);
            }
            else
            {
              Miss = Hbar[2];
              eMiss = eHbar[2];
            }
            Miss = 1.;

            Double_t xx = (Double_t)nmb / (Double_t)nclock;
            Double_t yy = ( nphoton - Miss*npion*2. ) / nev[id];
            Double_t eyy = sqrt( nphoton + pow(eMiss*npion*2.,2.) + pow(Miss*enpion*2.,2.) ) / nev[id];
            if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
            {
              gr[ig]->SetPoint(igp[ig], xx, yy);
              gr[ig]->SetPointError(igp[ig], 0., eyy);
              gr_run[ig]->SetPoint(igp[ig], runnumber, yy);
              gr_run[ig]->SetPointError(igp[ig], 0., eyy);
              igp[ig]++;
            }
          } // ic, is
      } // ipt, id

    delete f_ert;
    delete f_mb;
    irun++;
  }

  TFile *f_out = new TFile(Form("pileup/Pileup-photon-%d.root",process), "RECREATE");
  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t id=0; id<2; id++)
      for(Int_t ic=0; ic<2; ic++)
        for(Int_t is=0; is<2; is++)
        {
          Int_t ig = ipt*8+id*4+ic*2+is;
          mcw( ig, Form("proc%d-data%d-cond%d-pt%d-%d", process, id, ic*2+is, pTlow[id][ipt], pThigh[id][ipt]) );
          gr[ig]->Set(igp[ig]);
          gr_run[ig]->Set(igp[ig]);
          gr[ig]->Write();
          gr_run[ig]->Write();
        }
  f_out->Close();
}
