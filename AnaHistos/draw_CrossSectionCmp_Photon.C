#include "GlobalVars.h"
#include "ReadGraph.h"

void draw_CrossSectionCmp_Photon()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};

  double gx[4][npT] = {}, gy[4][npT] = {}, egy[4][npT] = {};
  for(int part=0; part<4; part++)
    ReadGraph<TGraphErrors>("data/CrossSection-photon.root", part, gx[part], gy[part], egy[part]);

  TGraph *gr_sasha = new TGraph("data/sasha-cross.txt");

  TGraphErrors *gr_parts[4];
  int igp_parts[4] = {};

  for(int part=0; part<4; part++)
  {
    gr_parts[part] = new TGraphErrors(npT);
    for(int ipt=2; ipt<npT; ipt++)
    {
      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double yy, eyy;

      int ip1 = Get_ipt(gx[part], xx);
      int ip2 = Get_ipt(gx[3], xx);

      double sasha_pT, sasha;
      gr_sasha->GetPoint(ipt-2, sasha_pT, sasha);

      if( TMath::Abs(sasha_pT - xx) > 0.2 )
      {
        cout << "Wrong pT matching!!!" << endl;
        return;
      }

      if(part < 3)
      {
        yy = ( gy[part][ip1] - gy[3][ip2] ) / gy[3][ip2];
        eyy = TMath::Abs(yy+1.) * sqrt( pow(egy[part][ip1]/gy[part][ip1],2) + pow(egy[3][ip2]/gy[3][ip2],2) );
      }
      else
      {
        yy = gy[3][ip2] / sasha;
        eyy = egy[3][ip2] / sasha;
      }

      if( eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_parts[part]->SetPoint(igp_parts[part], xx, yy);
        gr_parts[part]->SetPointError(igp_parts[part], 0., eyy);
        igp_parts[part]++;
      }
    }
    gr_parts[part]->Set(igp_parts[part]);
  }

  mc(0);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  mc(1);

  for(int part=0; part<4; part++)
  {
    mcd(part/3, 0);
    if(part < 3)
    {
      gr_parts[part]->SetTitle("Diff in parts;p_{T} [GeV];Diff;");
      aset(gr_parts[part], "","", 6.1,30., -0.5,0.5);
      leg0->AddEntry(gr_parts[part], Form("%s",pname[part]), "P");
    }
    else
    {
      gr_parts[part]->SetTitle("#gamma/#pi^{0};p_{T} [GeV];#gamma/#pi^{0};");
      aset(gr_parts[part], "","", 6.1,30.);
    }
    style(gr_parts[part], part+20, part+1);
    if(part%3 == 0)
      gr_parts[part]->Draw("APE");
    else
      gr_parts[part]->Draw("PE");
    if(part == 0)
      leg0->Draw();
  }

  c0->Print("plots/CrossSectionCmpParts-photon.pdf");
  c1->Print("plots/CrossSectionCmpCombined-photon.pdf");
}
