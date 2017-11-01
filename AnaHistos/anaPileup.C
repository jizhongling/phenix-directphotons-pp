#include "Pileup.h"
#include "BBCCounts.h"

void anaPileup(const Int_t process = 0)
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TGraphErrors *gr[npT*8];
  Int_t igp[npT*8] = {};
  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t id=0; id<2; id++)
      for(Int_t ic=0; ic<2; ic++)
        for(Int_t is=0; is<2; is++)
        {
          Int_t ig = ipt*8+id*4+ic*2+is;
          mc(ig, 5,4);
          gr[ig] = new TGraphErrors(20);
          gr[ig]->SetName(Form("gr_%d",ig));
        }

  const Int_t nThread = 20;
  Int_t thread = -1;
  Int_t irun = 0;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

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
    THnSparse *hn_pion[2];
    hn_pion[0] = (THnSparse*)f_ert->Get("hn_pion");
    hn_pion[1] = (THnSparse*)f_mb->Get("hn_pion");

    ULong64_t nclock = GetClockLive(runnumber);
    ULong64_t nmb = GetBBCNarrowLive(runnumber);
    Double_t nev[2];
    nev[0] = h_events_ert->GetBinContent( h_events_ert->GetXaxis()->FindBin("ert_c") );
    nev[1] = h_events_mb->GetBinContent( h_events_mb->GetXaxis()->FindBin("bbc_novtx_10cm") );

    TF1 *fn_fit = new TF1("fn_fit", "gaus+pol2(3)", 0.06, 0.25);
    TF1 *fn_bg = new TF1("fn_bg", "pol2", 0.06, 0.25);

    for(Int_t ipt=0; ipt<npT; ipt++)
      for(Int_t id=0; id<2; id++)
        for(Int_t ic=0; ic<2; ic++)
          for(Int_t is=0; is<2; is++)
          {
            Int_t ig = ipt*8+id*4+ic*2+is;
            mcd(ig, irun+1);

            hn_pion[id]->GetAxis(3)->SetRange(ic+1,ic+1);
            hn_pion[id]->GetAxis(0)->SetRange(secl[is],sech[is]);
            hn_pion[id]->GetAxis(1)->SetRange(pTlow[id][ipt], pThigh[id][ipt]);

            TH1 *h_minv = hn_pion[id]->Projection(2);
            h_minv->Rebin(10);
            Double_t max = h_minv->GetMaximum();
            if( max <= 0. )
            {
              delete h_minv;
              continue;
            }

            h_minv->SetTitle(Form("#%d",runnumber));
            aset(h_minv, "m_{inv} [GeV]","", 0.,0.3);

            Double_t par[10] = {max,0.140,0.010, 0.,0.,0.};
            fn_fit->SetParameters(par);
            h_minv->Fit(fn_fit, "RQ0");
            fn_fit->GetParameters(par);
            fn_bg->SetParameters(par+3);

            fn_fit->SetLineColor(kRed);
            fn_bg->SetLineColor(kGreen);
            h_minv->DrawCopy("EHIST");
            fn_fit->DrawCopy("SAME");
            fn_bg->DrawCopy("SAME");

            Double_t nsig = 0.;
            Double_t nbg = 0.;
            for(Int_t ib=12; ib<=16; ib++)
            {
              nsig += h_minv->GetBinContent(ib);
              Double_t bincenter = h_minv->GetXaxis()->GetBinCenter(ib);
              nbg += fn_bg->Eval(bincenter);
            }
            if( max < 5. )
              nbg = ( h_minv->Integral(6,10) + h_minv->Integral(19,23) ) / 2.;

            Double_t ensig = sqrt(nsig);
            Double_t rbg = nbg / nsig;
            Double_t erbg = sqrt(nbg) / nsig;
            Double_t npion = nsig * (1-rbg);
            Double_t enpion = sqrt( pow(ensig*(1-rbg),2.) + pow(nsig*erbg,2.) );

            Double_t xx = (Double_t)nmb / (Double_t)nclock;
            Double_t yy = npion / nev[id];
            Double_t eyy = enpion / nev[id];
            if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
            {
              gr[ig]->SetPoint(igp[ig], xx, yy);
              gr[ig]->SetPointError(igp[ig], 0., eyy);
              igp[ig]++;
            }
            delete h_minv;
          }

    delete fn_fit;
    delete fn_bg;
    delete f_ert;
    delete f_mb;
    irun++;
  }

  TFile *f_out = new TFile(Form("pileup/Pileup-%d.root",process), "RECREATE");
  for(Int_t ipt=0; ipt<npT; ipt++)
    for(Int_t id=0; id<2; id++)
      for(Int_t ic=0; ic<2; ic++)
        for(Int_t is=0; is<2; is++)
        {
          Int_t ig = ipt*8+id*4+ic*2+is;
          gr[ig]->Set(igp[ig]);
          gROOT->ProcessLine( Form("c%d->Print(\"pileup/Minv-proc%d-data%d-cond%d-pt%d-%d.pdf\");", ig, process, id, ic*2+is, pTlow[id][ipt], pThigh[id][ipt]) );
          gr[ig]->Write();
        }
  f_out->Close();
}
