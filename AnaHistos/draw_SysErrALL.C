#include "GlobalVars.h"
#include "QueryTree.h"
#include "IsoPhotonALL.h"

void draw_SysErrALL()
{
  const char *beam_list[3] = {"A_{L}^{Blue}", "A_{L}^{Yellow}", "A_{LL}"};

  QueryTree *qt_sys = new QueryTree("data/IsoPhotonALL-syserr.root", "RECREATE");
  QueryTree *qt_all = new QueryTree("data/IsoPhotonALL.root");

  mc();
  mcd();
  double estat, dummy;
  qt_all->Query(0, 2+ngr_photon*2, dummy, dummy, estat);
  for(int beam=0; beam<3; beam++)
  {
    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      double xpt, comb[2], ecomb[2];
      for(int isys=0; isys<2; isys++)
      {
        int igr = beam + ngr_photon*2 + ngr_photon*3*isys;
        qt_all->Query(ipt, igr, xpt, comb[isys], ecomb[isys]);
      }
      double ediff = beam==2 ? 2.*estat : 0.;
      double esys = sqrt(pow(comb[1]-comb[0],2) + pow(3.853e-4,2) + ediff*ediff)*1.066;
      qt_sys->Fill(ipt, beam, xpt, comb[0], esys);
    } // ipt

    int igr = beam + ngr_photon*2;
    TGraphErrors *gr_all = qt_all->Graph(igr);
    TGraphErrors *gr_sys = qt_sys->Graph(beam);
    gr_all->SetTitle( Form("#gamma^{dir} %s",beam_list[beam]) );
    aset(gr_all, "p_{T} [GeV]",beam_list[beam], 0.,20., -0.06,0.05);
    style(gr_all, 1, 1);
    style(gr_sys, 1, 1);
    gr_all->Draw("AP");
    gr_sys->Draw("[]");
    c0->Print(Form("plots/IsoPhotonALL-beam%d-syserr.pdf",beam));
    c0->Clear("D");
  } // beam

  qt_sys->Save();
}
