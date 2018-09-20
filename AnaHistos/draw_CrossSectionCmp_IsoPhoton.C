#include "GlobalVars.h"
#include "ReadGraph.h"

void draw_CrossSectionCmp_IsoPhoton()
{
  const double PI = TMath::Pi();
  const double jetphox_scale = 1./200.;  // combined 200 histograms
  const char *fname[3] = {"halfpt", "onept", "twopt"};
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};

  double gx[4][npT] = {}, gy[4][npT] = {}, egy[4][npT] = {};
  for(int part=0; part<4; part++)
    ReadGraph<TGraphErrors>("data/CrossSection-isophoton.root", part, gx[part], gy[part], egy[part]);

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  TGraphErrors *gr_parts[3];
  int igp_parts[3] = {};

  for(int part=0; part<3; part++)
  {
    gr_parts[part] = new TGraphErrors(npT);
    for(int ipt=12; ipt<npT; ipt++)
    {
      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;

      int ip1 = Get_ipt(gx[part], xx);
      int ip2 = Get_ipt(gx[3], xx);

      double yy = ( gy[part][ip1] - gy[3][ip2] ) / gy[3][ip2];
      double eyy = TMath::Abs(yy+1.) * sqrt( pow(egy[part][ip1]/gy[part][ip1],2) + pow(egy[3][ip2]/gy[3][ip2],2) );

      if( eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_parts[part]->SetPoint(igp_parts[part], xx, yy);
        gr_parts[part]->SetPointError(igp_parts[part], 0., eyy);
        igp_parts[part]++;
      }
    }
    gr_parts[part]->Set(igp_parts[part]);
  }

  TGraphErrors *gr_g2pi0;
  int ipg_g2pi0 = 0;

  gr_g2pi0 = new TGraphErrors(npT);
  for(int ipt=2; ipt<npT; ipt++)
  {
    double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
    double yy, eyy;

    int ip2 = Get_ipt(gx[3], xx);

    double sasha_pT, sasha;
    gr_sasha->GetPoint(ipt-2, sasha_pT, sasha);

    if( TMath::Abs(sasha_pT - xx) > 0.2 )
    {
      cout << "Wrong pT matching!!!" << endl;
      return;
    }

    yy = gy[3][ip2] / sasha;
    eyy = egy[3][ip2] / sasha;

    if( eyy > 0. && eyy < TMath::Infinity() )
    {
      gr_g2pi0->SetPoint(ipg_g2pi0, xx, yy);
      gr_g2pi0->SetPointError(ipg_g2pi0, 0., eyy);
      ipg_g2pi0++;
    }
  }
  gr_g2pi0->Set(ipg_g2pi0);

  TGraphErrors *gr_nlo[3];
  int igp_nlo[3] = {};

  for(int imu=0; imu<3; imu++)
  {
    gr_nlo[imu] = new TGraphErrors(npT);
    TFile *f_nlo = new TFile( Form("data/isoprompt-ct10-%s.root",fname[imu]) );
    TH1 *h_nlo = (TH1*)f_nlo->Get("hp41");
    h_nlo->Scale(jetphox_scale);

    for(int ipt=12; ipt<npT; ipt++)
    {
      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;

      double factor = 1. / (2*PI*xx*0.5);
      int bin_th = h_nlo->GetXaxis()->FindBin(xx);
      double nnlo = factor * h_nlo->GetBinContent(bin_th);
      double ennlo = factor * h_nlo->GetBinError(bin_th);

      int ip2 = Get_ipt(gx[3], xx);

      double yy = gy[3][ip2] / nnlo;
      double eyy = yy * sqrt( pow(egy[3][ip2]/gy[3][ip2],2) + pow(ennlo/nnlo,2) );

      if( eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_nlo[imu]->SetPoint(igp_nlo[imu], xx, yy);
        gr_nlo[imu]->SetPointError(igp_nlo[imu], 0., eyy);
        igp_nlo[imu]++;
      }
    }
    gr_nlo[imu]->Set(igp_nlo[imu]);

    delete h_nlo;
    delete f_nlo;
  }

  mc(0);
  mcd(0);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(int part=0; part<3; part++)
  {
    gr_parts[part]->SetTitle("Diff in parts;p_{T} [GeV];Diff;");
    aset(gr_parts[part], "","", 6.1,30., -0.5,0.5);
    leg0->AddEntry(gr_parts[part], Form("%s",pname[part]), "P");
    style(gr_parts[part], part+20, part+1);
    if(part == 0)
      gr_parts[part]->Draw("APE");
    else
      gr_parts[part]->Draw("PE");
  }
  leg0->Draw();

  mc(1);
  mcd(1);
  gr_g2pi0->SetTitle("Photon to pi0 ratio;p_{T} [GeV];#gamma/#pi^{0};");
  aset(gr_g2pi0, "","", 6.1,30.);
  style(gr_g2pi0, 20, 1);
  gr_g2pi0->Draw("APE");

  mc(2);
  mcd(2);
  legi(1, 0.2,0.8,0.9,0.9);
  leg1->SetNColumns(3);

  for(int imu=0; imu<3; imu++)
  {
    gr_nlo[imu]->SetTitle("data/theory;p_{T} [GeV];#frac{data}{theory}");
    aset(gr_nlo[imu], "","", 6.1,30., 0.5,2);
    leg1->AddEntry(gr_nlo[imu], Form("%s",fname[imu]), "L");
    style(gr_nlo[imu], imu+20, imu+1);
    if(imu == 0)
      gr_nlo[imu]->Draw("ALE");
    else
      gr_nlo[imu]->Draw("LE");
  }
  leg1->Draw();

  c0->Print("plots/CrossSectionCmpParts-isophoton.pdf");
  c1->Print("plots/CrossSectionCmp2pi0-isophoton.pdf");
  c2->Print("plots/CrossSectionCmp2Jetphox-isophoton.pdf");
}
