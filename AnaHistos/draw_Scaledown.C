#include "csc-runno-bbcsc.h"

void AddHistos(TH1 **h_out, THnSparse *hn_pion, int scaledown)
{
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  const int rangel[3] = {47, 112, 187};
  const int rangeh[3] = {97, 162, 227};
  const double rangew[3] = {-0.5, 1., -0.5};

  for(int part=0; part<2; part++)
  {
    hn_pion->GetAxis(0)->SetRange(secl[part],sech[part]);
    for(int ir=0; ir<3; ir++)
    {
      hn_pion->GetAxis(2)->SetRange(rangel[ir],rangeh[ir]);
      TH1 *h_tmp = hn_pion->Projection(1);
      TH1 *h_pion = (TH1*)h_tmp->Clone("h_pion");
      delete h_tmp;
      h_out[part]->Add(h_pion, rangew[ir] * (scaledown+1) );
    }
  }
}

void AddRun(TH1 **h_ertb, TH1 **h_ertc, TH1 **h_mb, int runno, int ertc_scaledown, int bbc_scaledown)
{
  TFile *f_ert = new TFile( Form("/phenix/plhf/zji/taxi/Run13pp510ERT/10678/data/DirectPhotonPP-%d.root",runno) );
  TFile *f_mb = new TFile( Form("/phenix/plhf/zji/taxi/Run13pp510MinBias/10701/data/DirectPhotonPP-%d.root",runno) );

  THnSparse *hn_pion_ert = (THnSparse*)f_ert->Get("hn_pion");
  THnSparse *hn_pion_mb = (THnSparse*)f_mb->Get("hn_pion");

  for(int ibit=2; ibit<=8; ibit+=2)
  {
    hn_pion_ert->GetAxis(5)->SetRange(ibit,ibit);
    AddHistos(h_ertb, hn_pion_ert, 0);
  }

  hn_pion_ert->GetAxis(5)->SetRange(5,8);
  AddHistos(h_ertc, hn_pion_ert, ertc_scaledown);

  hn_pion_mb->GetAxis(5)->SetRange(17,32);
  AddHistos(h_mb, hn_pion_mb, bbc_scaledown);

  return;
}

void draw_Scaledown()
{
  gROOT->ProcessLine(".L ReadGraph.C");

  const double pTbins[32] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,
    100.0 };

  double gx[30], TrigE[3][2][30], eTrigE[3][2][30];
  const int ispion = 1;
  for(int trig=0; trig<3; trig++)
    for(int part=0; part<2; part++)
      ReadGraphAsymmErrors("TriggerEfficiency.root", 9*ispion+3*trig+part, gx, (double*)TrigE[trig][part], (double*)eTrigE[trig][part]);

  const double TrigE_BBC = 0.91;
  const double eTrigE_BBC = 0.01;

  TH1 *h_sum_ertb[3];
  TH1 *h_sum_ertc[3];
  TH1 *h_sum_mb[3];
  TH1 *h_ertb[3][2];
  TH1 *h_ertc[3][2];
  TH1 *h_mb[3][2];

  for(int csc=0; csc<3; csc++)
  {
    h_sum_ertb[csc] = new TH1D(Form("h_ertb_csc%d",csc), "ERT_4x4b", 31, pTbins);
    h_sum_ertc[csc] = new TH1D(Form("h_ertc_csc%d",csc), "ERT_4x4c", 31, pTbins);
    h_sum_mb[csc] = new TH1D(Form("h_mb_csc%d",csc), "MinBias", 31, pTbins);
    for(int part=0; part<2; part++)
    {
      h_ertb[csc][part] = new TH1D(Form("h_ertb_csc%d_part%d",csc,part), "ERT_4x4b", 31, pTbins);
      h_ertc[csc][part] = new TH1D(Form("h_ertc_csc%d_part%d",csc,part), "ERT_4x4c", 31, pTbins);
      h_mb[csc][part] = new TH1D(Form("h_mb_csc%d_part%d",csc,part), "MinBias", 31, pTbins);
    }
  }

  for(int ir=0; ir<sizeof(csc0_run)/sizeof(csc0_run[0]); ir++)
    AddRun((TH1**)h_ertb[0], (TH1**)h_ertc[0], (TH1**)h_mb[0], csc0_run[ir], 0, csc0_bbc[ir]);
  for(int ir=0; ir<sizeof(csc1_run)/sizeof(csc1_run[0]); ir++)
    AddRun((TH1**)h_ertb[1], (TH1**)h_ertc[1], (TH1**)h_mb[1], csc1_run[ir], 1, csc1_bbc[ir]);
  for(int ir=0; ir<sizeof(csc2_run)/sizeof(csc2_run[0]); ir++)
    AddRun((TH1**)h_ertb[2], (TH1**)h_ertc[2], (TH1**)h_mb[2], csc2_run[ir], 2, csc2_bbc[ir]);

  //for(int ir=0; ir<2; ir++)
  //{
  //  AddRun((TH1**)h_ertb[0], (TH1**)h_ertc[0], (TH1**)h_mb[0], csc0_run[ir], 0, csc0_bbc[ir]);
  //  AddRun((TH1**)h_ertb[1], (TH1**)h_ertc[1], (TH1**)h_mb[1], csc1_run[ir], 1, csc1_bbc[ir]);
  //  AddRun((TH1**)h_ertb[2], (TH1**)h_ertc[2], (TH1**)h_mb[2], csc2_run[ir], 2, csc2_bbc[ir]);
  //}

  for(int csc=0; csc<3; csc++)
    for(int ipt=1; ipt<30; ipt++)
    {
      double npion_ertb = 0.;
      double npion_ertc = 0.;
      double npion_mb = 0.;
      for(int part=0; part<2; part++)
      {
        npion_ertb += h_ertb[csc][part]->GetBinContent(ipt+1) / TrigE[1][part][ipt];
        npion_ertc += h_ertc[csc][part]->GetBinContent(ipt+1) / TrigE[2][part][ipt];
        npion_mb += h_mb[csc][part]->GetBinContent(ipt+1) / TrigE_BBC;
      }
      h_sum_ertb[csc]->SetBinContent(ipt+1, npion_ertb);
      h_sum_ertc[csc]->SetBinContent(ipt+1, npion_ertc);
      h_sum_mb[csc]->SetBinContent(ipt+1, npion_mb);
    }

  TCanvas *c = new TCanvas("c", "", 1200, 1200);
  gStyle->SetOptStat(0);
  c->Divide(2,2);
  int ipad = 1;

  for(int csc=0; csc<3; csc++)
  {
    c->cd(ipad++);
    gPad->SetLogy();
    h_sum_ertc[csc]->SetTitle( Form("4x4c scaledown %d",csc) );
    h_sum_ertc[csc]->GetXaxis()->SetRange(1., 30.);
    h_sum_ertb[csc]->SetLineColor(1);
    h_sum_ertc[csc]->SetLineColor(2);
    h_sum_mb[csc]->SetLineColor(3);
    h_sum_ertc[csc]->Draw();
    h_sum_ertb[csc]->Draw("SAME");
    h_sum_mb[csc]->Draw("SAME");
  }

  c->cd(ipad++);
  TLegend *leg = new TLegend(0.2, 0.2, 0.8, 0.8);
  leg->AddEntry(h_sum_ertb[0], "4x4b", "L");
  leg->AddEntry(h_sum_ertc[0], "4x4c", "L");
  leg->AddEntry(h_sum_mb[0], "MB", "L");
  leg->Draw();

  c->Print("plots/Scaledown.pdf");
}
