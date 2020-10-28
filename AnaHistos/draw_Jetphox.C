#include "GlobalVars.h"

void draw_Jetphox()
{
  const char *jetphox_fname[5] = {"MMM", "MLL", "MLH", "MHL", "MHH"};

  mc();
  mcd();
  legi(0, 0.3,0.25,0.7,0.35);
  leg0->SetNColumns(2);

  double sigma_onept[npT], esigma_onept[npT];
  for(int imu=0; imu<5; imu++)
  {
    TGraphErrors *gr_ratio = new TGraphErrors(npT);
    TFile *f_nlo = new TFile( Form("data/isoprompt-x%d-ct14-%s.root",imu==0?2000:200,jetphox_fname[imu]) );
    TH1 *h_nlo = (TH1*)f_nlo->Get("hp41");

    int igp = 0;
    for(int ipt=12; ipt<npT; ipt++)
    {
      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      int bin_th = h_nlo->GetXaxis()->FindBin(xpt);
      double sigma_nlo = h_nlo->GetBinContent(bin_th);
      double esigma_nlo = h_nlo->GetBinError(bin_th);

      if(imu == 0)
      {
        sigma_onept[ipt] = sigma_nlo / 10.;
        esigma_onept[ipt] = esigma_nlo / 10.;
      }
      else
      {
        double ratio = sigma_nlo/sigma_onept[ipt];
        double eratio = ratio*sqrt(pow(esigma_onept[ipt]/sigma_onept[ipt],2) + pow(esigma_nlo/sigma_nlo,2));
        gr_ratio->SetPoint(igp, xpt, ratio);
        gr_ratio->SetPointError(igp, 0., eratio);
        igp++;
      }
    } // ipt

    if(imu == 0)
      continue;

    gr_ratio->Set(igp);
    style(gr_ratio, 1, imu);
    aset(gr_ratio, "p_{T} [GeV/c]","Ratio", 6.,30., 0.8,1.2);
    gr_ratio->Draw(imu==1?"AL":"L");
    leg0->AddEntry(gr_ratio, Form("%s",jetphox_fname[imu]), "L");
  } // imu
  leg0->Draw();

  c0->Print("plots/Jetphox.pdf");
}
