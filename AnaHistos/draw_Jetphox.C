#include "GlobalVars.h"
#include "QueryTree.h"

void draw_Jetphox()
{
  const double jetphox_scale = 1./400.;  // combined 400 histograms
  const char *jetphox_fname[3] = {"halfpt", "onept", "twopt"};
  const char *jetphox_setnum[1] = {"isoprompt"};
  const char *jetphox_setden[1] = {"incprompt"};

  QueryTree *qt_ratio = new QueryTree("data/JetphoxRatio.root", "RECREATE");

  mc(0);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(int iset=0; iset<1; iset++)
  {
    mcd(0, iset+1);
    for(int imu=0; imu<3; imu++)
    {
      TFile *f_num = new TFile( Form("data/%s-x400-ct14-%s.root",jetphox_setnum[iset],jetphox_fname[imu]) );
      TFile *f_den = new TFile( Form("data/%s-x400-ct14-%s.root",jetphox_setden[iset],jetphox_fname[imu]) );
      TH1 *h_num = (TH1*)f_num->Get("hp41");
      TH1 *h_den = (TH1*)f_den->Get("hp41");
      h_num->Scale(jetphox_scale);
      h_den->Scale(jetphox_scale);

      for(int ipt=10; ipt<npT; ipt++)
      {
        double xpt = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
        int bin_th = h_num->GetXaxis()->FindBin(xpt);
        double ratio = h_num->GetBinContent(bin_th) / h_den->GetBinContent(bin_th);
        double eratio = ratio * sqrt( pow(h_num->GetBinError(bin_th)/h_num->GetBinContent(bin_th),2) + pow(h_den->GetBinError(bin_th)/h_den->GetBinContent(bin_th),2) );
        if( TMath::Finite(ratio+eratio) )
          qt_ratio->Fill(ipt, imu, xpt, ratio, eratio);
      }

      TGraphErrors *gr_ratio = qt_ratio->Graph(imu);
      gr_ratio->SetTitle( Form("%s/%s;p_{T} [GeV/c];#frac{%s}{%s}",jetphox_setnum[iset],jetphox_setden[iset],jetphox_setnum[iset],jetphox_setden[iset]) );
      if(iset < 2)
      {
        aset(gr_ratio, "","", 5.1,30., 0.7,1.1);
      }
      else if(iset == 2)
      {
        aset(gr_ratio, "","", 5.1,30., 0.,0.3);
      }
      else if(iset == 3)
      {
        aset(gr_ratio, "","", 5.1,30., 0.,0.6);
      }
      style(gr_ratio, imu+20, imu+1);
      if(imu == 0)
        gr_ratio->Draw("ALE");
      else
        gr_ratio->Draw("LE");
      if(iset == 0)
        leg0->AddEntry(gr_ratio, Form("%s",jetphox_fname[imu]), "L");

      delete f_num;
      delete f_den;
    }
  }
  leg0->Draw();

  //c0->Print("plots/JetphoxRatio.pdf");
  qt_ratio->Save();
}
