#include "BBCCounts.h"

void anaPileup(const Int_t process = 0)
{
  TGraphErrors *gr[4];
  for(Int_t ig=0; ig<4; ig++)
  {
    mc(ig, 5,4);
    gr[ig] = new TGraphErrors(20);
    gr[ig]->SetName(Form("gr_%d",ig));
  }

  const Int_t nThread = 20;
  Int_t thread = -1;
  Int_t irun = 0;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");

  ReadClockCounts();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    const Int_t secl[2] = {1, 7};
    const Int_t sech[2] = {6, 8};

    THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
    TAxis *axis = hn_pion->GetAxis(2);
    Int_t bin112 = axis->FindBin(0.112);
    Int_t bin162 = axis->FindBin(0.162);

    TF1 *fn_init = new TF1("fn_init", "gaus", 0., 0.3);
    TF1 *fn_sig = new TF1("fn_sig", "gaus(0)+pol2(3)", 0., 0.3);
    TF1 *fn_bg = new TF1("fn_bg", "pol2", 0., 0.3);

    ULong64_t nclock = GetClockLive(runnumber);
    ULong64_t nmb = GetBBCNarrowLive(runnumber);
    ULong64_t scaledown = GetBBCNovtxScaledown(runnumber) + 1; 

    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<2; is++)
      {
        mcd(ic*2+is, irun+1);

        hn_pion->GetAxis(3)->SetRange(ic+1,ic+1);
        hn_pion->GetAxis(0)->SetRange(secl[is],sech[is]);
        hn_pion->GetAxis(1)->SetRange(5,30);
        TH1 *h_minv = hn_pion->Projection(2);
        aset(h_minv, "m_{inv} [GeV]","", 0.,0.3);

        Double_t par[10];
        h_minv->Fit(fn_init, "Q0", "", 0.112, 0.162);
        fn_sig->SetParameters( fn_init->GetParameters() );
        h_minv->Fit(fn_sig, "Q0", "", 0.047, 0.227);
        fn_sig->GetParameters(par);
        fn_bg->SetParameters(par[3], par[4], par[5]);

        fn_bg->SetLineColor(kGreen);
        h_minv->DrawCopy();
        fn_bg->Draw("SAME");

        Double_t nsig = 0.;
        Double_t nbg = 0.;
        for(Int_t ib=bin112; ib<bin162; ib++)
        {
          nsig += h_minv->GetBinContent(ib);
          Double_t bincenter = axis->GetBinCenter(ib);
          nbg += fn_bg->Eval(bincenter);
        }
        Double_t npion = nsig - nbg;
        Double_t enpion = sqrt(1./nsig + 1./nbg);

        Double_t xx = (Double_t)nmb/(Double_t)nclock;
        Double_t yy = npion * (Double_t)scaledown / (Double_t)nmb;
        Double_t eyy = enpion * (Double_t)scaledown / (Double_t)nmb;
        gr[ic*2+is]->SetPoint(irun, xx, yy);
        gr[ic*2+is]->SetPointError(irun, 0., eyy);
        delete h_minv;
      }

    delete f;
    irun++;
  }

  TFile *f_out = new TFile(Form("pileup/Mine-%d.root",process), "RECREATE");
  for(Int_t ig=0; ig<4; ig++)
  {
    gROOT->ProcessLine( Form("c%d->Print(\"pileup/Mine-proc%d-cond%d.pdf\");", ig, process, ig) );
    gr[ig]->Write();
  }
  f_out->Close();
}
