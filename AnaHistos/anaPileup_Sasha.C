#include "BBCCounts.h"

void anaPileup_Sasha(const Int_t process = 0)
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
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  ReadClockCounts();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510MinBias/12148/data/Pi0PP-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("hevtype");
    TH1 *h_minv[2][2];  // h_minv[ic][is]
    Double_t npion[2][2];  // npion[ic][is]
    Double_t enpion[2][2];  // enpion[ic][is]
    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<2; is++)
      {
        h_minv[ic][is] = new TH1F(Form("h_minv_%d",ic*2+is), Form("#%d",runnumber), 1000,0.,1.);
        npion[ic][is] = 0.;
        enpion[ic][is] = 0.;
      }

    const char *cname[2] = {"p", "tp"};
    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<3; is++)
        for(Int_t ip=4; ip<20; ip++)
        {
          TH1 *mchist = (TH1*)f->Get(Form("mc_s%d_bcc0_pt_%03d_%s",is,5*ip,cname[ic]));
          h_minv[ic][is/2]->Add(mchist);
        }

    TF1 *fn_fit = new TF1("fn_fit", "gaus+pol2(3)", 0.06, 0.25);
    TF1 *fn_bg = new TF1("fn_bg", "pol2", 0.06, 0.25);

    ULong64_t nclock = GetClockLive(runnumber);
    ULong64_t nmb = GetBBCNarrowLive(runnumber);
    Double_t nev = h_events->GetBinContent(1);

    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<2; is++)
      {
        mcd(ic*2+is, irun+1);
        h_minv[ic][is]->Rebin(10);
        h_minv[ic][is]->Scale(0.5);
        aset(h_minv[ic][is], "m_{inv} [GeV]","", 0.,0.3);

        Double_t par[10] = {h_minv[ic][is]->GetMaximum(),0.140,0.010, 0.,0.,0.};
        fn_fit->SetParameters(par);
        h_minv[ic][is]->Fit(fn_fit, "RQ0");
        fn_fit->GetParameters(par);
        fn_bg->SetParameters(par+3);

        fn_fit->SetLineColor(kRed);
        fn_bg->SetLineColor(kGreen);
        h_minv[ic][is]->DrawCopy("EHIST");
        fn_fit->DrawCopy("SAME");
        fn_bg->DrawCopy("SAME");

        Double_t nsig = 0.;
        Double_t nbg = 0.;
        for(Int_t ib=12; ib<=16; ib++)
        {
          nsig += h_minv[ic][is]->GetBinContent(ib);
          Double_t bincenter = h_minv[ic][is]->GetXaxis()->GetBinCenter(ib);
          nbg += fn_bg->Eval(bincenter);
        }
        npion[ic][is] = nsig - nbg;
        enpion[ic][is] = sqrt(nsig + nbg);
      }

    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<2; is++)
      {
        Double_t xx = (Double_t)nmb/(Double_t)nclock;
        Double_t yy = npion[ic][is] / nev;
        Double_t eyy = enpion[ic][is] / nev;
        if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
        {
          gr[ic*2+is]->SetPoint(irun, xx, yy);
          gr[ic*2+is]->SetPointError(irun, 0., eyy);
        }
        delete h_minv[ic][is];
      }
    delete f;
    irun++;
  }

  TFile *f_out = new TFile(Form("pileup/Sasha-%d.root",process), "RECREATE");
  for(Int_t ig=0; ig<4; ig++)
  {
    gROOT->ProcessLine( Form("c%d->Print(\"pileup/Sasha-proc%d-cond%d.pdf\");", ig, process, ig) );
    gr[ig]->Write();
  }
  f_out->Close();
}
