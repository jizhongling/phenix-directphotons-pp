#include "GlobalVars.h"
#include "ReadGraph.h"

void draw_CrossSectionCmp()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};

  Double_t gx[4][npT] = {}, gy[4][npT] = {}, egy[4][npT] = {};
  for(Int_t part=0; part<4; part++)
    ReadGraph<TGraphErrors>("CrossSection-pion.root", part, gx[part], gy[part], egy[part]);

  TGraph *gr_sasha = new TGraph("sasha-cross.txt");
  const Double_t scale_sasha = 1e9;

  TGraphErrors *gr_ratio[4];
  Int_t igp[4] = {};

  for(Int_t part=0; part<4; part++)
  {
    gr_ratio[part] = new TGraphErrors(npT);

    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t yy, eyy;

      Int_t ip1 = Get_ipt(gx[part], xx);
      Int_t ip2 = Get_ipt(gx[3], xx);

      if(part < 3)
      {
        yy = ( gy[part][ip1] - gy[3][ip2] ) / gy[3][ip2];
        eyy = TMath::Abs(yy+1.) * sqrt( pow(egy[part][ip1]/gy[part][ip1],2.) + pow(egy[3][ip2]/gy[3][ip2],2.) );
      }
      else
      {
        Double_t Sasha = gr_sasha->Eval(xx) * scale_sasha;
        yy = ( gy[3][ip2] - Sasha ) / Sasha;
        eyy = egy[3][ip2] / Sasha;
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

  mc(0, 2,1);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(Int_t part=0; part<4; part++)
  {
    mcd(0, part/3+1);
    if(part <3)
    {
      gr_ratio[part]->SetTitle("Diff in parts;p_{T} [GeV];Diff;");
      aset(gr_ratio[part], "","", 6.,20., -0.2,0.2);
      leg0->AddEntry(gr_ratio[part], Form("%s",pname[part]), "P");
    }
    else
    {
      gr_ratio[part]->SetTitle("Diff with Sasha;p_{T} [GeV];Diff;");
      aset(gr_ratio[part], "","", 6.,20., -0.3,0.);
    }
    style(gr_ratio[part], part+20, part+1);
    if(part%3 == 0)
      gr_ratio[part]->Draw("APE");
    else
      gr_ratio[part]->Draw("PE");
  }
  mcd(0, 1);
  leg0->Draw();

  c0->Print("CrossSectionCmp.pdf");
}
