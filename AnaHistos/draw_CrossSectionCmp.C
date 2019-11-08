#include "GlobalVars.h"
#include "QueryTree.h"

void draw_CrossSectionCmp(const int nameid)
{
  TString name;
  switch(nameid):
  {
    case 0:
      name = "pion";
      break;
    case 1:
      name = "photon";
      break;
    case 2:
      name = "isophoton";
      break;
    default:
      cout << "Wrong nameid" << endl;
      return;
  }

  const double PI = TMath::Pi();
  const double DeltaEta = 0.5;
  const double jetphox_scale = 1./400.;  // combined 400 histograms
  const char *jetphox_fname[3] = {"halfpt", "onept", "twopt"};
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};

  TString fname = "data/CrossSection-";
  fname += name + ".root";
  QueryTree *qt_cross = new QueryTree(fname);

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  TGraphErrors *gr_parts[4];
  int igp_parts[4] = {};
  for(int part=0; part<4; part++)
    gr_parts[part] = new TGraphErrors(npT);

  for(int ipt=2; ipt<npT; ipt++)
  {
    double xpt, Combine, eCombine;
    if( !qt_cross->Query(ipt, 3, xpt, Combine, eCombine) )
      continue;

    for(int part=0; part<3; part++)
    {
      double Cross, eCross;
      qt_cross->Query(ipt, part, xpt, Cross, eCross);
      double yy = (Cross - Combine) / Combine;
      double eyy = TMath::Abs(yy+1.) * sqrt( pow(eCross/Cross,2) + pow(eCombine/Combine,2) );
      if( TMath::Finite(yy+eyy) )
      {
        gr_parts[part]->SetPoint(igp_parts[part], xpt, yy);
        gr_parts[part]->SetPointError(igp_parts[part], 0., eyy);
        igp_parts[part]++;
      }
    }

    double sasha_pT, sasha;
    gr_sasha->GetPoint(ipt-2, sasha_pT, sasha);
    if( TMath::Abs(sasha_pT - xpt) > 0.2 )
    {
      cout << "ipt " << ipt << ": wrong pT matching!!!" << endl;
      //return;
    }

    double yy = Combine / sasha;
    double eyy = eCombine / sasha;
    if( name.EqualTo("pion") )
      yy -= 1.;
    if( TMath::Finite(yy+eyy) )
    {
      gr_parts[part]->SetPoint(igp_parts[part], xpt, yy);
      gr_parts[part]->SetPointError(igp_parts[part], 0., eyy);
      igp_parts[part]++;
    }
  }

  mc(0);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  mc(1);

  for(int part=0; part<4; part++)
  {
    gr_parts[part]->Set(igp_parts[part]);

    mcd(part/3, 0);
    if(part < 3)
    {
      gr_parts[part]->SetTitle("Diff in parts;p_{T} [GeV];Diff;");
      aset(gr_parts[part], "","", 0.,30., -0.5,0.5);
      leg0->AddEntry(gr_parts[part], Form("%s",pname[part]), "P");
    }
    else
    {
      if( name.EqualTo("pion") )
      {
        gr_parts[part]->SetTitle("Diff with PRD 93, 011501;p_{T} [GeV];Diff;");
        aset(gr_parts[part], "","", 6.1,30., -0.1,0.1);
      }
      else
      {
        gr_parts[part]->SetTitle("#gamma/#pi^{0};p_{T} [GeV];#gamma/#pi^{0};");
        aset(gr_parts[part], "","", 6.1,30.);
      }
    }
    style(gr_parts[part], part+20, part+1);

    if(part%3 == 0)
      gr_parts[part]->Draw("APE");
    else
      gr_parts[part]->Draw("PE");
    if(part == 0)
      leg0->Draw();
  }

  fname = "plots/CrossSectionCmpParts-";
  fname += name + ".pdf";
  c0->Print(fname);
  fname = "plots/CrossSectionCmpCombined-";
  fname += name + ".pdf";
  c1->Print(fname);

  if( !name.EqualTo("isophoton") )
    return;

  mc(2);
  mcd(2);
  legi(1, 0.2,0.8,0.9,0.9);
  leg1->SetNColumns(3);

  TLine *line = new TLine();
  line->SetLineColor(kRed);
  line->SetLineWidth(5);
  line->SetLineStyle(2);

  for(int imu=0; imu<3; imu++)
  {
    TGraphErrors *gr_nlo = new TGraphErrors(npT);
    TFile *f_nlo = new TFile( Form("data/isoprompt-x400-ct14-%s.root",jetphox_fname[imu]) );
    TH1 *h_nlo = (TH1*)f_nlo->Get("hp41");
    h_nlo->Scale(jetphox_scale);

    int igp = 0;
    for(int ipt=12; ipt<npT; ipt++)
    {
      double xpt, Combine, eCombine;
      if( !qt_cross->Query(ipt, 3, xpt, Combine, eCombine) )
        continue;

      double factor = 1. / (2*PI*xpt*DeltaEta);
      int bin_th = h_nlo->GetXaxis()->FindBin(xpt);
      double sigma_nlo = factor * h_nlo->GetBinContent(bin_th);
      double esigma_nlo = factor * h_nlo->GetBinError(bin_th);

      double yy = Combine / sigma_nlo;
      double eyy = yy * sqrt( pow(eCombine/Combine,2) + pow(esigma_nlo/sigma_nlo,2) );
      if( TMath::Finite(yy+eyy) )
      {
        gr_nlo->SetPoint(igp, xpt, yy);
        gr_nlo->SetPointError(igp, 0., eyy);
        igp++;
      }
    }

    gr_nlo->Set(igp);
    gr_nlo->SetTitle("data/theory;p_{T} [GeV];#frac{data}{theory}");
    aset(gr_nlo, "","", 6.1,30., 0.5,2);
    leg1->AddEntry(gr_nlo, Form("%s",jetphox_fname[imu]), "L");
    style(gr_nlo, imu+20, imu+1);
    if(imu == 0)
      gr_nlo->Draw("ALE");
    else
      gr_nlo->Draw("LE");
    line->DrawLine(6.1, 1., 30., 1.);

    delete h_nlo;
    delete f_nlo;
  }
  leg1->Draw();

  c2->Print("plots/CrossSectionCmp2Jetphox-isophoton.pdf");
}
