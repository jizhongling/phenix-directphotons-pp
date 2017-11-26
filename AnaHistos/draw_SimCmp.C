#include "DivideFunctions.h"

void GenerateAcceptance(TFile *fsig, TFile *ftot, TObjArray *Glist, Int_t isim)
{
  const Int_t npT = 31;
  const Double_t vpT[npT] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0 };

  TF1 *cross = new TF1("cross", "x*(1/(1+exp((x-[5])/[6]))*[0]/pow(1+x/[1],[2])+(1-1/(1+exp((x-[5])/[6])))*[3]/pow(x,[4]))", 0, 30);
  cross->SetParameters(2.02819e+04, 4.59173e-01, 7.51170e+00, 1.52867e+01, 7.22708e+00, 2.15396e+01, 3.65471e+00);

  TH1::SetDefaultSumw2();

  TH1 *h_sig[3];
  for(Int_t part=0; part<3; part++)
    h_sig[part] = new TH1D(Form("h_sig_%d",part), "#pi^{0} signal count; p_{T} [GeV/c];", npT-1,vpT);
  TH1 *h_tot;

  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};
  char *pname[3] = {"PbScW", "PbScE", "PbGlE"};  // must be non-const if used directly in TLegend

  THnSparse *hn_total = (THnSparse*)fsig->Get("hn_total");
  //hn_total->GetAxis(3)->SetRange(0,80);
  for(Int_t part=0; part<3; part++)
  {
    hn_total->GetAxis(2)->SetRange(secl[part],sech[part]);
    for(Int_t ipt=0; ipt<30; ipt++)
    {
      hn_total->GetAxis(0)->SetRange(ipt+1, ipt+1);
      TH1 *h_tmp = hn_total->Projection(0);
      Double_t npion = h_tmp->GetEntries();
      h_sig[part]->SetBinContent(ipt+1, npion);
      h_sig[part]->SetBinError( ipt+1, sqrt(npion) * cross->Eval(vpT[ipt]) );
      delete h_tmp;
    }
    h_tot = (TH1*)ftot->Get("h_pion");
  }

  mc(1);
  legi(1, 0.2,0.7,0.5,0.9);

  for(Int_t part=0; part<3; part++)
  {
    TGraphErrors *gr = DivideHisto(h_sig[part], h_tot);
    //TGraphAsymmErrors *gr = new TGraphAsymmErrors(h_sig[part], h_tot, "n");
    Glist->Add(gr);
    mcd(1);
    aset(gr, "p_{T} [GeV]", "acceptance", 0.,30., 0.,0.16);
    style(gr, 20+part, 1+part);
    gr->SetName(Form("gr_%d",3*isim+part));
    gr->SetTitle("#pi^{0} acceptance");
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    leg1->AddEntry(gr, pname[part], "P");
    leg1->Draw();
  }

  c1->Print(Form("plots/SimCmp-%d.pdf",isim));
  delete c1;

  for(Int_t part=0; part<3; part++)
    delete h_sig[part];
  delete h_tot;
  delete hn_total;

  return;
}

void DrawCmp(TObjArray *Glist)
{
  Double_t gx[30];
  Double_t gy[6][30];
  Double_t egy[6][30];
  for(Int_t isim=0; isim<2; isim++)
    for(Int_t part=0; part<3; part++)
      ReadGraph((TGraph*)Glist->At(3*isim+part), gx, (Double_t*)gy[3*isim+part], (Double_t*)egy[3*isim+part]);

  mc();
  legi(0, 0.2,0.7,0.5,0.9);

  TGraphErrors *gr_ratio[3];
  Double_t rgy[3][30] = {};
  Double_t ergy[3][30] = {};

  for(Int_t part=0; part<3; part++)
  {
    for(Int_t ipt=0; ipt<30; ipt++)
    {
      rgy[part][ipt] = ( gy[3+part][ipt] - gy[part][ipt] ) / gy[part][ipt];
      ergy[part][ipt] = (gy[3+part][ipt]/gy[part][ipt]) * sqrt( pow(egy[3+part][ipt]/gy[3+part][ipt],2.) + pow(egy[part][ipt]/gy[part][ipt],2.) );
    }
    gr_ratio[part] =  new TGraphErrors(30, gx, rgy[part], 0, ergy[part]);
  }

  const char *pname[3] = {"PbScW", "PbScE", "PbGlE"};
  for(Int_t part=0; part<3; part++)
  {
    mcd();
    aset(gr_ratio[part]);
    style(gr_ratio[part], 20+part, 1+part);
    gr_ratio[part]->SetTitle("Ratio;p_{T} [GeV];ratio;");
    gr_ratio[part]->GetXaxis()->SetRangeUser(0., 30.);
    gr_ratio[part]->GetYaxis()->SetRangeUser(-0.3, 0.3);
    if(part==0)
      gr_ratio[part]->Draw("APE");
    else
      gr_ratio[part]->Draw("PE");
    leg0->AddEntry(gr_ratio[part], Form("%s",pname[part]), "P");
  }
  leg0->Draw();

  c0->Print("plots/SimCmp.pdf");
  delete c0;

  return;
}

void draw_SimCmp()
{
  gROOT->ProcessLine(".L ReadGraph.C");

  TFile *fsig = new TFile("data/MissingRatio-histo.root");
  TFile *ftot = new TFile("data/AnaPHPythia-histo.root");
  TObjArray *Glist = new TObjArray();

  GenerateAcceptance(ftot, ftot, Glist, 0);
  GenerateAcceptance(fsig, ftot, Glist, 1);

  DrawCmp(Glist);

  TFile *fout = new TFile("data/SimCmp.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
