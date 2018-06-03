#include "GlobalVars.h"
#include "ReadGraph.h"

void draw_CrossSectionCmp_Photon()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};

  Double_t gx[4][npT] = {}, gy[4][npT] = {}, egy[4][npT] = {};
  for(Int_t part=0; part<4; part++)
    ReadGraph<TGraphErrors>("data/CrossSection-photon.root", part, gx[part], gy[part], egy[part]);

  TGraphErrors *gr_ratio[3];
  Int_t igp[3] = {};

  for(Int_t part=0; part<3; part++)
  {
    gr_ratio[part] = new TGraphErrors(npT);

    for(Int_t ipt=2; ipt<npT; ipt++)
    {
      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t yy, eyy;

      Int_t ip1 = Get_ipt(gx[part], xx);
      Int_t ip2 = Get_ipt(gx[3], xx);

      yy = ( gy[part][ip1] - gy[3][ip2] ) / gy[3][ip2];
      eyy = TMath::Abs(yy+1.) * sqrt( pow(egy[part][ip1]/gy[part][ip1],2.) + pow(egy[3][ip2]/gy[3][ip2],2.) );

      if( eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_ratio[part]->SetPoint(igp[part], xx, yy);
        gr_ratio[part]->SetPointError(igp[part], 0., eyy);
        igp[part]++;
      }
    }

    gr_ratio[part]->Set(igp[part]);
  }

  mc();
  mcd();
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  leg0->Draw();

  for(Int_t part=0; part<3; part++)
  {
    gr_ratio[part]->SetTitle("Diff in parts;p_{T} [GeV];Diff;");
    aset(gr_ratio[part], "","", 6.1,30., -0.5,0.5);
    leg0->AddEntry(gr_ratio[part], Form("%s",pname[part]), "P");
    style(gr_ratio[part], part+20, part+1);
    if(part == 0)
      gr_ratio[part]->Draw("APE");
    else
      gr_ratio[part]->Draw("PE");
  }

  c0->Print("plots/CrossSectionCmpParts-photon.pdf");
}

#include "GlobalVars.h"
#include "ReadGraph.h"

void draw_CrossSectionCmp_Photon()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};

  Double_t gx[4][npT] = {}, gy[4][npT] = {}, egy[4][npT] = {};
  for(Int_t part=0; part<4; part++)
    ReadGraph<TGraphErrors>("data/CrossSection-photon.root", part, gx[part], gy[part], egy[part]);

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  TGraphErrors *gr_ratio[4];
  Int_t igp[4] = {};

  for(Int_t part=0; part<4; part++)
  {
    gr_ratio[part] = new TGraphErrors(npT);

    for(Int_t ipt=2; ipt<npT; ipt++)
    {
      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t yy, eyy;

      Int_t ip1 = Get_ipt(gx[part], xx);
      Int_t ip2 = Get_ipt(gx[3], xx);

      Double_t sasha_pT, sasha;
      gr_sasha->GetPoint(ipt-2, sasha_pT, sasha);

      if( TMath::Abs(sasha_pT - xx) > 0.2 )
      {
        cout << "Wrong pT matching!!!" << endl;
        return;
      }

      if(part < 3)
      {
        yy = ( gy[part][ip1] - gy[3][ip2] ) / gy[3][ip2];
        eyy = TMath::Abs(yy+1.) * sqrt( pow(egy[part][ip1]/gy[part][ip1],2.) + pow(egy[3][ip2]/gy[3][ip2],2.) );
      }
      else
      {
        yy = gy[3][ip2]/ sasha;
        eyy = egy[3][ip2] / sasha;
      }

      if( eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_ratio[part]->SetPoint(igp[part], xx, yy);
        gr_ratio[part]->SetPointError(igp[part], 0., eyy);
        igp[part]++;
      }
    }

    gr_ratio[part]->Set(igp[part]);
  }

  mc(0);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  leg0->Draw();
  mc(1);

  for(Int_t part=0; part<4; part++)
  {
    mcd(part/3, 0);
    if(part < 3)
    {
      gr_ratio[part]->SetTitle("Diff in parts;p_{T} [GeV];Diff;");
      aset(gr_ratio[part], "","", 6.1,30., -0.5,0.5);
      leg0->AddEntry(gr_ratio[part], Form("%s",pname[part]), "P");
    }
    else
    {
      gr_ratio[part]->SetTitle("#gamma/#pi^{0};p_{T} [GeV];#gamma/#pi^{0};");
      aset(gr_ratio[part], "","", 6.1,30., 0.,0.65);
    }
    style(gr_ratio[part], part+20, part+1);
    if(part%3 == 0)
      gr_ratio[part]->Draw("APE");
    else
      gr_ratio[part]->Draw("PE");
  }

  c0->Print("plots/CrossSectionCmpParts-photon.pdf");
  c1->Print("plots/CrossSectionCmpCombined-photon.pdf");
}
