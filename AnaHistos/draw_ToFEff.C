#include "GlobalVars.h"

Bool_t FitMinv(THnSparse *hn_pion, Int_t ipt, Double_t &npion, Double_t &enpion)
{
  TF1 *fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", 0.06, 0.25);
  TF1 *fn_bg = new TF1("fn_bg", "pol2", 0.06, 0.25);

  hn_pion->GetAxis(1)->SetRange(ipt+1,ipt+1);
  TH1 *h_minv = hn_pion->Projection(2); 
  h_minv->Rebin(10);
  Double_t max = h_minv->GetMaximum();
  if( max <= 0. )
  {
    delete h_minv;
    return kFALSE;
  }

  h_minv->SetTitle(Form("pT: %4.2f-%4.2f",pTbin[ipt],pTbin[ipt+1]));
  aset(h_minv, "m_{inv} [GeV]","", 0.,0.3);

  Double_t par[10] = {max,0.140,0.010, 0.,0.,0.};
  for(Int_t ifit=0; ifit<5; ifit++)
  {
    fn_fit->SetParameters(par);
    h_minv->Fit(fn_fit, "RQ0");
    fn_fit->GetParameters(par);
  }
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

  Int_t ndf = fn_fit->GetNDF();
  Double_t prob = fn_fit->GetProb();
  if( ndf < 10 || prob < 0.1 ) 
    nbg = ( h_minv->Integral(6,10) + h_minv->Integral(19,23) ) / 2.;

  Double_t ensig = sqrt(nsig);
  Double_t rbg = nbg / nsig;
  Double_t erbg = sqrt(nbg) / nsig;

  npion = nsig * (1-rbg);
  enpion = sqrt( pow(ensig*(1-rbg),2.) + pow(nsig*erbg,2.) );

  delete h_minv;
  return kTRUE;
}

void draw_ToFEff()
{
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};
  const char *name[2] = {"PbSc", "PbGl"};

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonNode-histo.root");

  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

  TGraphErrors *gr_tof[2];
  Int_t igr[2] = {};
  for(Int_t part=0; part<2; part++)
    gr_tof[part] =  new TGraphErrors(npT);

  mc(0, 5,6);
  mc(1, 5,6);
  mc(2, 2,1);

  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      hn_pion->GetAxis(0)->SetRange(secl[part],sech[part]);

      mcd(0, ipt+1);
      Double_t nt, ent;
      hn_pion->GetAxis(3)->SetRange(1,1);
      FitMinv(hn_pion, ipt, nt, ent);

      mcd(1, ipt+1);
      Double_t np, enp;
      hn_pion->GetAxis(3)->SetRange(2,2);
      FitMinv(hn_pion, ipt, np, enp);

      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t yy = np / nt;
      Double_t eyy = yy * sqrt( pow(ent/nt,2.) + pow(enp/np,2.) );
      if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_tof[part]->SetPoint(igr[part], xx, yy);
        gr_tof[part]->SetPointError(igr[part], 0., eyy);
        igr[part]++;
      }
    }

  for(Int_t part=0; part<2; part++)
  {
    mcd(2, part+1);
    gr_tof[part]->SetTitle( Form("ToF efficeincy for %s", name[part]) );
    aset(gr_tof[part], "p_{T} [GeV]", "Eff", 0.,30., 0.,1.1);
    style(gr_tof[part], 24, kRed);
    gr_tof[part]->Draw("APE");
    gr_tof[part]->Fit("pol0", "Q","", 9.,30.);

    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr_tof[part]->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetY1NDC(0.6);
    st->SetX2NDC(1.0);
    st->SetY2NDC(0.8);
  }

  c0->Print("ToFEff-PbSc.pdf");
  c1->Print("ToFEff-PbGl.pdf");
  c2->Print("ToFEff.pdf");
}
