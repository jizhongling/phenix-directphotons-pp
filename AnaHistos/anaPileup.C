#include "BBCCounts.h"

void anaPileup(const Int_t process = 0)
{
  TGraphErrors *gr[4];
  TGraphErrors *gr_run[4];
  for(Int_t ig=0; ig<4; ig++)
  {
    mc(ig, 5,4);
    gr[ig] = new TGraphErrors(20);
    gr[ig]->SetName(Form("gr_%d",ig));
    gr_run[ig] = new TGraphErrors(20);
    gr_run[ig]->SetName(Form("gr_run_%d",ig));
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

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    const Int_t secl[2] = {1, 7};
    const Int_t sech[2] = {6, 8};

    TH1 *h_events = (TH1*)f->Get("h_events");
    THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

    TF1 *fn_fit = new TF1("fn_fit", "gaus+pol2(3)", 0.06, 0.25);
    TF1 *fn_bg = new TF1("fn_bg", "pol2", 0.06, 0.25);

    ULong64_t nclock = GetClockLive(runnumber);
    ULong64_t nmb = GetBBCNarrowLive(runnumber);
    Double_t nev = h_events->GetBinContent(1);

    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<2; is++)
      {
        mcd(ic*2+is, irun+1);

        hn_pion->GetAxis(3)->SetRange(ic+1,ic+1);
        hn_pion->GetAxis(0)->SetRange(secl[is],sech[is]);
        hn_pion->GetAxis(1)->SetRange(5,20);
        TH1 *h_minv = hn_pion->Projection(2);
        h_minv->Rebin(10);
        h_minv->SetTitle(Form("#%d",runnumber));
        aset(h_minv, "m_{inv} [GeV]","", 0.,0.3);

        Double_t par[10] = {h_minv->GetMaximum(),0.140,0.010, 0.,0.,0.};
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
        Double_t npion = nsig - nbg;
        Double_t enpion = sqrt(nsig + nbg);

        Double_t xx = (Double_t)nmb / (Double_t)nclock;
        Double_t yy = npion / nev;
        Double_t eyy = enpion / nev;
        if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
        {
          gr[ic*2+is]->SetPoint(irun, xx, yy);
          gr[ic*2+is]->SetPointError(irun, 0., eyy);
          gr_run[ic*2+is]->SetPoint(irun, runnumber, yy);
          gr_run[ic*2+is]->SetPointError(irun, 0., eyy);
        }
        delete h_minv;
      }

    delete f;
    irun++;
  }

  TFile *f_out = new TFile(Form("pileup/Mine-%d.root",process), "RECREATE");
  TFile *f_out = new TFile("Pileup.root", "RECREATE");
  for(Int_t ig=0; ig<4; ig++)
  {
    gROOT->ProcessLine( Form("c%d->Print(\"pileup/Mine-proc%d-cond%d.pdf\");", ig, process, ig) );
    gr[ig]->Write();
    //gr_run[ig]->Write();
  }
  f_out->Close();
}
