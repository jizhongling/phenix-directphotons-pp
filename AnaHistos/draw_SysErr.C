#include "GlobalVars.h"
#include "QueryTree.h"
#include "Chi2Fit.h"

void draw_SysErr()
{
  for(int iso=0; iso<2; iso++)
  {
    char *type = iso ? "iso" : "";
    QueryTree *qt_sys = new QueryTree(Form("data/CrossSection-%sphoton-syserr.root",type), "RECREATE");
    QueryTree *qt_cross = new QueryTree(Form("data/CrossSection-%sphoton.root",type));

    for(int ipt=0; ipt<npT; ipt++)
    {
      double xpt, xsec, exsec, rsys[2], ersys[2], rfit[3], erfit[3], ren[3], eren[3];
      qt_cross->Query(ipt, 3, xpt, xsec, exsec);
      for(int part=0; part<3; part++)
      {
        qt_cross->Query(ipt, 1+part+3*1, xpt, rfit[part], erfit[part]);
        qt_cross->Query(ipt, 1+part+3*2, xpt, ren[part], eren[part]);
      }
      rsys[0] = TMath::MaxElement(3, rfit);
      Chi2Fit(3, ren, eren, rsys[1], ersys[1]);
      double sys = xsec*sqrt(rsys[0]*rsys[0] + rsys[1]*rsys[1]);
      qt_sys->Fill(ipt, 0, xpt, xsec, sys);
    }

    mc();
    mcd();
    gPad->SetLogy();
    TGraphErrors *gr_cross = qt_cross->Graph(3);
    TGraphErrors *gr_sys = qt_sys->Graph(0);
    aset(gr_cross, "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{-3}]", 6.1,30., 1e-1, 2e3);
    style(gr_cross, 1, 1);
    style(gr_sys, 1, 1);
    gr_cross->Draw("AP");
    gr_sys->Draw("[]");
    c0->Print(Form("plots/CrossSection-%sphoton-syserr.pdf",type));
    qt_sys->Save();
  }
}
