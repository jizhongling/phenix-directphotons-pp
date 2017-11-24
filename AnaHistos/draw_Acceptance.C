#include "GlobalVars.h"
#include "FitMinv.h"
#include "GetEfficiency.h"

void draw_Acceptance()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  TF1 *cross = new TF1("cross", "x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))", 0, 30);
  cross->SetParameters(2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00);

  TGraphAsymmErrors *gr[3];
  Int_t igp[3] = {};
  for(Int_t part=0; part<3; part++)
  {
    gr[part] = new TGraphAsymmErrors(npT);
    gr[part]->SetName(Form("gr_%d",part));
    mc(part, 6,5);
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-warn-histo.root");

  TH1 *h_total = (TH1*)f->Get("h_pion");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  TAxis *axis_pt0 = hn_pion->GetAxis(0);
  TAxis *axis_pt = hn_pion->GetAxis(1);
  TAxis *axis_minv = hn_pion->GetAxis(2);
  TAxis *axis_sec = hn_pion->GetAxis(3);
  TAxis *axis_peak = hn_pion->GetAxis(4);

  for(Int_t part=0; part<3; part++)
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t ww = cross->Eval(xx);

      TH1 *h_tt = (TH1*)h_total->Clone();
      h_tt->Scale(1./ww);
      Double_t nt = h_tt->GetBinContent(ipt+1);
      Double_t ent = sqrt(nt);
      delete h_tt;

      Double_t np, enp;
      mcd(part, ipt+1);
      //axis_peak->SetRange(3,3);  // Require two peaks
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(ipt+1,ipt+1);
      TH1 *h_minv = hn_pion->Projection(2);
      h_minv->Rebin(10);
      h_minv->Scale(1./ww);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV",pTbin[ipt],pTbin[ipt+1]) );
      FitMinv(h_minv, np, enp);
      // Do not use fit result, use all invariant mass instead
      np = h_minv->Integral(0,-1);
      enp = sqrt(np);
      delete h_minv;

      Double_t yy, eyyl, eyyh;
      if( !GetEfficiency(nt,np, yy,eyyl,eyyh) )
      {
        eyyl = yy * sqrt( pow(ent/nt,2.) + pow(enp/np,2.) );
        eyyh = 0.;
      }
      if( yy >= 0. && eyyl >= 0. && eyyl < TMath::Infinity() )
      {
        gr[part]->SetPoint(igp[part], xx, yy);
        gr[part]->SetPointError(igp[part], 0.,0., eyyl,eyyh);
        igp[part]++;
      }
    }


  mc(3);
  mcd(3);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(Int_t part=0; part<3; part++)
  {
    gr[part]->Set(igp[part]);
    gr[part]->SetTitle("#pi^{0} acceptance");
    aset(gr[part], "p_{T} [GeV]","acceptance", 0.,30., 0.,0.12);
    style(gr[part], part+20, part+1);
    if(part==0)
      gr[part]->Draw("AP");
    else
      gr[part]->Draw("P");
    leg0->AddEntry(gr[part], pname[part], "P");
    TGraph *gr_sasha =  new TGraph( Form("sasha-acc-part%d.txt",part) );
    gr_sasha->Draw("C");
  }

  leg0->Draw();
  c3->Print("Acceptance.pdf");

  TFile *f_out = new TFile("data/Acceptance.root", "RECREATE");
  for(Int_t part=0; part<3; part++)
  {
    mcw( part, Form("part%d",part) );
    gr[part]->Write();
  }
  f_out->Close();
}
