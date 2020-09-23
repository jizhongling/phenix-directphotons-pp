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
      double yy = Cross/Combine;
      double eyy = yy*sqrt(pow(eCross/Cross,2) + pow(eCombine/Combine,2));
      if( TMath::Finite(yy+eyy) )
      {
        gr_parts[part]->SetPoint(igp_parts[part], xpt, yy-1);
        gr_parts[part]->SetPointError(igp_parts[part], 0., eyy);
        igp_parts[part]++;
      }
    } // part

    double sasha_pT, sasha;
    gr_sasha->GetPoint(ipt-2, sasha_pT, sasha);
    if( TMath::Abs(sasha_pT - xpt) > 0.2 )
    {
      cout << "ipt " << ipt << ": wrong pT matching!!!" << endl;
      //return;
    }

    double yy = Combine/sasha;
    double eyy = eCombine/sasha;
    if( TMath::Finite(yy+eyy) )
    {
      gr_parts[3]->SetPoint(igp_parts[3], xpt, nameid>0?yy:yy-1);
      gr_parts[3]->SetPointError(igp_parts[3], 0., eyy);
      igp_parts[3]++;
    }
  } // ipt

  mc(0);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  mc(1);

  for(int part=0; part<4; part++)
  {
    gr_parts[part]->Set(igp_parts[part]);

    if(part < 3)
    {
      mcd(0);
      gr_parts[part]->SetTitle("Diff in parts;p_{T} [GeV/c];Diff;");
      aset(gr_parts[part], "","", 5.9,30.1, -0.5,0.5);
      leg0->AddEntry(gr_parts[part], Form("%s",pname[part]), "P");
    }
    else
    {
      mcd(1);
      if( name.EqualTo("pion") )
      {
        gr_parts[part]->SetTitle("Diff with PRD 93, 011501;p_{T} [GeV/c];Diff;");
        aset(gr_parts[part], "","", 5.9,30.1, -0.05,0.05);
      }
      else
      {
        gr_parts[part]->SetTitle("#gamma/#pi^{0};p_{T} [GeV/c];#gamma/#pi^{0};");
        aset(gr_parts[part], "","", 5.9,30.1, 0.,0.5);
      }
    }
    style(gr_parts[part], part+20, part+1);

    char *opt = part%3==0 ? "APE" : "PE";
    gr_parts[part]->Draw(opt);
    if(part == 0)
      leg0->Draw();
  }

  fname = "plots/CrossSectionCmpParts-";
  fname += name + ".pdf";
  c0->Print(fname);
  fname = "plots/CrossSectionCmpCombined-";
  fname += name + ".pdf";
  c1->Print(fname);
}
