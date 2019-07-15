#include "GlobalVars.h"
#include "QueryTree.h"

void draw_Iso2Incl_Data()
{
  const int sector = 3;  // PbSc west: 0; PbSc east: 1; PbGl: 2; Combined: 3

  QueryTree *qt_pion = new QueryTree("data/CrossSection-pion.root");
  QueryTree *qt_photon = new QueryTree("data/CrossSection-photon.root");
  QueryTree *qt_isophoton = new QueryTree("data/CrossSection-isophoton.root");

  TGraphErrors *gr[2];
  int igp[2] = {};
  for(int iph=0; iph<2; iph++)
    gr[iph] = new TGraphErrors(npT);

  for(int ipt=0; ipt<npT; ipt++)
  {
    double xpt, incl, eincl, iso, eiso;

    if( qt_pion->Query(ipt, sector, xpt, incl, eincl) &&
        qt_pion->Query(ipt, sector+4, xpt, iso, eiso) )
    {
      double yy = iso / incl;
      double eyy = yy * sqrt( pow(eincl/incl,2) + pow(eiso/iso,2) );
      if( TMath::Finite(yy+eyy) )
      {
        gr[0]->SetPoint(igp[0], xpt, yy);
        gr[0]->SetPointError(igp[0], 0., eyy);
        igp[0]++;
      }
    }

    if( qt_photon->Query(ipt, sector, xpt, incl, eincl) &&
        qt_isophoton->Query(ipt, sector, xpt, iso, eiso) )
    {
      double yy = iso / incl;
      double eyy = yy * sqrt( pow(eincl/incl,2) + pow(eiso/iso,2) );
      if( TMath::Finite(yy+eyy) )
      {
        gr[1]->SetPoint(igp[1], xpt, yy);
        gr[1]->SetPointError(igp[1], 0., eyy);
        igp[1]++;
      }
    }
  }

  mc();
  mcd();
  legi(0, 0.2,0.7,0.4,0.9);
  for(int iph=0; iph<2; iph++)
  {
    gr[iph]->Set(igp[iph]);
    aset(gr[iph], "p_{T} [GeV]", "Iso/All", 6.,30., 0.,1.4);
    style(gr[iph], iph+20, iph+1);
    if(iph==0)
      gr[iph]->Draw("AP");
    else
      gr[iph]->Draw("P");
  }
  leg0->AddEntry(gr[0], "#pi^{0}", "P");
  leg0->AddEntry(gr[1], "#gamma_{dir}", "P");
  leg0->Draw();
  c0->Print("plots/Iso2Incl-data.pdf");
}
