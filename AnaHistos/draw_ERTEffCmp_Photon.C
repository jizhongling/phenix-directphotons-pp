#include "QueryTree.h"

void draw_ERTEffCmp_Photon()
{
  const char *pname[2] = {"PbSc", "PbGl"};

  QueryTree *qt_ert = new QueryTree("data/ERTEff-photon.root");
  QueryTree *qt_ert_mb = new QueryTree("data/ERTEff-photon-MB.root");

  mc();
  mcd();
  legi(0, 0.7,0.2,0.9,0.4);

  for(int part=0; part<2; part++)
  {
    TGraphAsymmErrors *gr = qt_ert->GraphAsymm(part);
    TGraphAsymmErrors *gr_mb = qt_ert_mb->GraphAsymm(part);
    gr->SetTitle("ERT_4x4c trigger efficeincy for photon");
    aset(gr, "p_{T} [GeV]","Eff", 1.,30., 0.,1.1);
    style(gr, part+20, part+1);
    style(gr_mb, part+24, part+1);
    if(part==0)
    {
      gr->Draw("APE");
      gr_mb->Draw("PE");
    }
    else
    {
      gr->Draw("PE");
      gr_mb->Draw("PE");
    }
    leg0->AddEntry(gr, pname[part], "LPE");
  }
  leg0->Draw();

  c0->Print("plots/ERTEffCmp-photon-MB.pdf");
}
