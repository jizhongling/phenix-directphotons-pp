#include "GlobalVars.h"

void draw_Jetphox()
{
  const double PI = TMath::Pi();
  const double jetphox_scale = 1./200.;  // combined 200 histograms
  const char *jetphox_fname[3] = {"halfpt", "onept", "twopt"};
  const char *jetphox_setnum[4] = {"isoprompt", "isoprompt-200gev", "onef", "onef-incl"};
  const char *jetphox_setden[4] = {"inclprompt", "inclprompt-200gev", "isoprompt", "inclprompt"};

  mc(0, 2,2);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(int iset=0; iset<4; iset++)
  {
    mcd(0, iset+1);
    for(int imu=0; imu<3; imu++)
    {
      TGraphErrors *gr_ratio = new TGraphErrors(npT);
      TFile *f_num = new TFile( Form("data/%s-ct10-%s.root",jetphox_setnum[iset],jetphox_fname[imu]) );
      TFile *f_den = new TFile( Form("data/%s-ct10-%s.root",jetphox_setden[iset],jetphox_fname[imu]) );
      TH1 *h_num = (TH1*)f_num->Get("hp41");
      TH1 *h_den = (TH1*)f_den->Get("hp41");
      h_num->Scale(jetphox_scale);
      h_den->Scale(jetphox_scale);

      int igp = 0;
      for(int ipt=12; ipt<npT; ipt++)
      {
        double xpt = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
        double factor = 1. / (2*PI*xpt*0.5);
        int bin_th = h_num->GetXaxis()->FindBin(xpt);
        double ratio = h_num->GetBinContent(bin_th) / h_den->GetBinContent(bin_th);
        double eratio = ratio * sqrt( pow(h_num->GetBinError(bin_th)/h_num->GetBinContent(bin_th),2) + pow(h_den->GetBinError(bin_th)/h_den->GetBinContent(bin_th),2) );

        if( TMath::Finite(ratio+eratio) )
        {
          gr_ratio->SetPoint(igp, xpt, ratio);
          gr_ratio->SetPointError(igp, 0., eratio);
          igp++;
        }
      }

      gr_ratio->Set(igp);
      gr_ratio->SetTitle( Form("%s/%s;p_{T} [GeV];#frac{%s}{%s}",jetphox_setnum[iset],jetphox_setden[iset],jetphox_setnum[iset],jetphox_setden[iset]) );
      if(iset < 2)
      {
        aset(gr_ratio, "","", 6.1,30., 0.6,1.5);
      }
      else if(iset == 2)
      {
        aset(gr_ratio, "","", 6.1,30., 0.,0.3);
      }
      else if(iset == 3)
      {
        aset(gr_ratio, "","", 6.1,30., 0.,0.6);
        leg0->AddEntry(gr_ratio, Form("%s",jetphox_fname[imu]), "L");
      }
      style(gr_ratio, imu+20, imu+1);
      if(imu == 0)
        gr_ratio->Draw("ALE");
      else
        gr_ratio->Draw("LE");

      delete h_num;
      delete h_den;
      delete f_num;
      delete f_den;
    }
  }
  leg0->Draw();

  c0->Print("plots/JetphoxRatio.pdf");
}
