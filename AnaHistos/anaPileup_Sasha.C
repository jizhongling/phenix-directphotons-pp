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
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");

  ReadClockCounts();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510MinBias/11955/data/Pi0PP-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *mh[2][2];  // mh[ic][is]
    Double_t npion[2][2];  // npion[ic][is]
    Double_t enpion[2][2];  // enpion[ic][is]
    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<2; is++)
      {
        mh[ic][is] = new TH1F(Form("mh_%d",ic*2+is), Form("#%d",runnumber), 1000,0.,1.);
        npion[ic][is] = 0.;
        enpion[ic][is] = 0.;
      }

    const char *cname[2] = {"p", "tp"};
    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<3; is++)
        for(Int_t ip=4; ip<40; ip++)
        {
          TH1 *mchist = (TH1*)f->Get(Form("mc_s%d_bcc0_pt_%03d_%s",is,5*ip,cname[ic]));
          mh[ic][is/2]->Add(mchist);
        }

    TAxis *axis = mh[0][0]->GetXaxis();
    Int_t bin112 = axis->FindBin(0.112);
    Int_t bin162 = axis->FindBin(0.162);

    TF1 *fn_init = new TF1("fn_init", "gaus", 0., 0.3);
    TF1 *fn_sig = new TF1("fn_sig", "gaus(0)+pol2(3)", 0., 0.3);
    TF1 *fn_bg = new TF1("fn_bg", "pol2", 0., 0.3);

    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<2; is++)
      {
        mcd(ic*2+is, irun+1);
        aset(mh[ic][is], "m_{inv} [GeV]","", 0.,0.3);

        Double_t par[10];
        mh[ic][is]->Fit(fn_init, "Q0", "", 0.112, 0.162);
        fn_sig->SetParameters( fn_init->GetParameters() );
        mh[ic][is]->Fit(fn_sig, "Q0", "", 0.047, 0.227);
        fn_sig->GetParameters(par);
        fn_bg->SetParameters(par[3], par[4], par[5]);

        fn_bg->SetLineColor(kGreen);
        mh[ic][is]->DrawCopy();
        fn_bg->Draw("SAME");

        Double_t nsig = 0.;
        Double_t nbg = 0.;
        for(Int_t ib=bin112; ib<bin162; ib++)
        {
          nsig += mh[ic][is]->GetBinContent(ib);
          Double_t bincenter = axis->GetBinCenter(ib);
          nbg += fn_bg->Eval(bincenter);
        }

        npion[ic][is] = nsig - nbg;
        enpion[ic][is] = sqrt(1./nsig + 1./nbg);
      }

    ULong64_t nclock = GetClockLive(runnumber);
    ULong64_t nmb = GetBBCNarrowLive(runnumber);
    ULong64_t scaledown = GetBBCNovtxScaledown(runnumber) + 1; 

    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<2; is++)
      {
        Double_t xx = (Double_t)nmb/(Double_t)nclock;
        Double_t yy = npion[ic][is] * (Double_t)scaledown / (Double_t)nmb;
        Double_t eyy = enpion[ic][is] * (Double_t)scaledown / (Double_t)nmb;
        gr[ic*2+is]->SetPoint(irun, xx, yy);
        gr[ic*2+is]->SetPointError(irun, 0., eyy);
        delete mh[ic][is];
      }
    delete f;
    irun++;
  }

  TFile *f_out = new TFile(Form("pileup/Pileup_Sasha-%d.root",process), "RECREATE");
  for(Int_t ig=0; ig<4; ig++)
  {
    gROOT->ProcessLine( Form("c%d->Print(\"pileup/Fit-proc%d-cond%d.pdf\");", ig, process, ig) );
    gr[ig]->Write();
  }
  f_out->Close();
}
