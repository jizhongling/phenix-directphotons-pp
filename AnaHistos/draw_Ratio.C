#include "GlobalVars.h"
#include "QueryTree.h"
#include "DivideFunctions.h"
#include "IsoPhotonALL.h"

void draw_Ratio()
{
  QueryTree *qt_num = new QueryTree("data/IsoPhotonALL-first3rd.root");
  QueryTree *qt_den = new QueryTree("data/IsoPhotonALL-last3rd.root");

  int beam = 2;
  int igr = beam + ngr_photon*2;
  TGraphErrors *gr_num = qt_num->Graph(igr);
  TGraphErrors *gr_den = qt_den->Graph(igr);
  TGraphErrors *gr_ratio = DivideGraph(gr_num, gr_den);

  mc();
  mcd();
  gPad->SetLogy();
  aset(gr_ratio, "p_{T} [GeV]","Ratio", 0.,30., 1e-2,1e2);
  style(gr_ratio, 20, 1);
  gr_ratio->Draw("AP");
  c0->Print("plots/ALLRatio.pdf");
}
