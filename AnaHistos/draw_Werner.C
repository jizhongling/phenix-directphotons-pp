#include "GlobalVars.h"
#include "QueryTree.h"
#include "Chi2Fit.h"

void draw_Werner()
{
  const double PI = TMath::Pi();
  const double DeltaEta = 0.5;

  const double jetphox_scale = 1./2000.;  // combined 2000 histograms
  const char *jetphox_fname[3] = {"onept", "halfpt", "twopt"};

  TGraphErrors *gr_cross, *gr_cross_sys;
  QueryTree *qt_cross = new QueryTree("data/CrossSection-isophoton.root");
  QueryTree *qt_sys = new QueryTree("data/CrossSection-syserr.root");

  TGraph *gr_werner = new TGraph("data/werner-cross.txt");

  TCanvas *c0 = new TCanvas("c0", "c0", 600,800);
  TPad *pad1 = new TPad("pad1", "pad1", 0.,0.35,1.,1., -1,0);
  TPad *pad2 = new TPad("pad2", "pad2", 0.,0.,1.,0.35, -1,0);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.);
  gPad->SetLogy();
  legi(0, 0.22,0.03,0.47,0.20);
  leg0->SetTextSize(0.035);

  pad2->cd();
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.25);

  TLine *line = new TLine();
  line->SetLineWidth(2);

  double nlo[npT][3], enlo[npT][3];
  for(int imu=0; imu<3; imu++)
  {
    TFile *f_nlo = new TFile( Form("data/isoprompt-x2000-ct14-%s.root",jetphox_fname[imu]) );
    TH1 *h_nlo = (TH1*)f_nlo->Get("hp41");
    h_nlo->Scale(jetphox_scale);
    for(int ipt=12; ipt<npT; ipt++)
    {
      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      int bin_th = h_nlo->GetXaxis()->FindBin(xpt);
      nlo[ipt][imu] = h_nlo->GetBinContent(bin_th);
      enlo[ipt][imu] = h_nlo->GetBinError(bin_th);
    }
    delete f_nlo;
  } // imu

  TGraphErrors *gr_central;
  for(int imu=0; imu<3; imu++)
  {
    TGraphErrors *gr_nlo = new TGraphErrors(npT);
    TGraphErrors *gr_ratio = new TGraphErrors(npT);
    TGraphErrors *gr_ratio_sys = new TGraphErrors(npT);

    int igp = 0;
    for(int ipt=12; ipt<npT; ipt++)
    {
      double xpt, Combine, eCombine, sysCombine;
      if( !qt_cross->Query(ipt, 3, xpt, Combine, eCombine) ||
          !qt_sys->Query(ipt, 1, xpt, Combine, sysCombine) )
        continue;

      double factor = 1. / (2*PI*xpt*DeltaEta);
      double sigma_nlo, esigma_nlo;
      sigma_nlo = nlo[ipt][imu] * factor;
      esigma_nlo = enlo[ipt][imu] * factor;
      gr_nlo->SetPoint(igp, xpt, sigma_nlo);
      gr_nlo->SetPointError(igp, 0., esigma_nlo);

      double ratio, eratio, sysratio;
      if(imu == 0)
      {
        ratio = Combine/sigma_nlo;
        eratio = ratio*sqrt(pow(eCombine/Combine,2) + pow(esigma_nlo/sigma_nlo,2));
        sysratio = sysCombine/sigma_nlo; 
        gr_ratio_sys->SetPoint(igp, xpt, ratio-1);
        gr_ratio_sys->SetPointError(igp, 0., sysratio);
      }
      else
      {
        ratio = sigma_nlo/nlo[ipt][0]/factor;
        eratio = ratio*sqrt(pow(enlo[ipt][0]/nlo[ipt][0],2) + pow(esigma_nlo/sigma_nlo,2));
      }
      gr_ratio->SetPoint(igp, xpt, ratio-1);
      gr_ratio->SetPointError(igp, 0., eratio);
      igp++;
    } // ipt

    gr_nlo->Set(igp);
    gr_ratio->Set(igp);
    gr_ratio_sys->Set(igp);

    style(gr_nlo, imu+1, imu+1, 2);
    style(gr_ratio, imu==0?20:imu+1, imu+1, 2);
    gr_ratio->SetMarkerSize(0.8);
    if(imu == 0)
      gr_central == gr_nlo;
    else
      leg0->AddEntry(gr_nlo, jetphox_fname[imu], "L");
    if(imu == 1)
      leg0->AddEntry(gr_central, jetphox_fname[0], "L");

    if(imu == 0)
    {
      pad1->cd();
      gr_cross = qt_cross->Graph(3);
      gr_cross_sys = qt_sys->Graph(1);
      gr_cross->SetTitle("");
      aset(gr_cross, "p_{T} [GeV/c]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{3}]", 4.9,30.1, 0.5e-1, 2e3);
      style(gr_cross, 20, 1, 2);
      style(gr_cross_sys, 1, 1, 2);
      style(gr_werner, 20, 3, 3);
      gr_cross->SetMarkerSize(0.8);
      gr_cross->Draw("AP");
      gr_cross_sys->Draw("[]");
      gr_nlo->Draw("LX");
      gr_werner->Draw("LX");
      leg0->AddEntry(gr_werner, "Werner", "L");
      leg0->Draw();

      pad2->cd();
      gr_ratio->SetTitle(";p_{T} [GeV/c];#frac{Data-Theory}{Theory}");
      aset(gr_ratio, "","", 4.9,30.1, -0.25,0.55, 1.,0.6,0.1,0.12);
      style(gr_ratio_sys, 1, imu+1, 2);
      gr_ratio->GetXaxis()->SetLabelSize(0.09);
      gr_ratio->GetYaxis()->SetLabelSize(0.09);
      gr_ratio->GetYaxis()->SetLabelOffset(0.01);
      gr_ratio->GetXaxis()->SetTickSize(0.08);
      gr_ratio->Draw("APE");
      gr_ratio_sys->Draw("[]");
      line->DrawLine(6., 0., 30., 0.);
    }
    else
    {
      pad1->cd();
      gr_nlo->Draw("LX");

      pad2->cd();
      gr_ratio->Draw("LX");
    }
  } // imu

  c0->Print("plots/Werner-isophoton.pdf");
}
