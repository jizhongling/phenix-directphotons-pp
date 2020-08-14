#include "GlobalVars.h"
#include "QueryTree.h"
#include "IsoPhotonALL.h"

void draw_SysErrALL()
{
  const char *beam_list[3] = {"A_{L}^{Blue}", "A_{L}^{Yellow}", "A_{LL}"};

  QueryTree *qt_sys = new QueryTree("data/IsoPhotonALL-syserr-tightcut.root", "RECREATE");
  QueryTree *qt_all = new QueryTree("data/IsoPhotonALL-tightcut.root");

  TBox *box = new TBox();
  box->SetLineColor(2);
  box->SetFillStyle(0);

  mc();
  mcd();
  gPad->SetGridy();
  for(int beam=0; beam<3; beam++)
  {
    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      double xpt, comb[3], ecomb[3];
      for(int isys=0; isys<3; isys++)
      {
        int igr = beam + ngr_photon*2 + ngr_photon*3*isys;
        qt_all->Query(ipt, igr, xpt, comb[isys], ecomb[isys]);
      }
      double esys = sqrt(pow(3.853e-4,2) + pow(comb[0]*0.066/(beam<2?2:1),2) + pow(comb[2]-comb[0],2) + pow(comb[1]-comb[0],2));
      qt_sys->Fill(ipt, beam, xpt, comb[0], esys);
    } // ipt

    int igr = beam + ngr_photon*2;
    TGraphErrors *gr_all = qt_all->Graph(igr);
    TGraphErrors *gr_sys = qt_sys->Graph(beam);
    gr_all->SetTitle( Form("#gamma^{dir} %s",beam_list[beam]) );
    aset(gr_all, "p_{T} [GeV]",beam_list[beam], 0.,20., -0.06,0.05);
    style(gr_all, 1, 1);
    style(gr_sys, 1, 1);
    gr_all->GetYaxis()->SetNdivisions(510);
    gr_sys->SetLineWidth(2);
    gr_all->Draw("AP");
    //gr_sys->Draw("[]");
    for(int i=0; i<gr_sys->GetN(); i++)
    {
      double xx, yy;
      gr_sys->GetPoint(i, xx, yy);
      double eyy = gr_sys->GetErrorY(i);
      box->DrawBox(xx-0.2,yy-eyy,xx+0.2,yy+eyy);
    }
    c0->Print(Form("plots/IsoPhotonALL-beam%d-syserr.pdf",beam));
    c0->Clear("D");
  } // beam

  qt_sys->Save();
}
