#include "GlobalVars.h"
#include "QueryTree.h"
#include "IsoPhotonALL.h"

void draw_SysErrALL()
{
  const char *beam_list[3] = {"A_{L}^{Blue}", "A_{L}^{Yellow}", "A_{LL}"};

  QueryTree *qt_sys = new QueryTree("data/IsoPhotonALL-syserr.root", "RECREATE");
  QueryTree *qt_all = new QueryTree("data/IsoPhotonALL.root");

  TGraph *gr_dssv = new TGraph("data/dssv-all.txt", "%lg %lg");

  TBox *box = new TBox();
  box->SetLineColor(2);
  box->SetFillStyle(0);

  const int nge = 25;
  Double_t xge[nge], yge[nge], eyge[nge];
  for(int i=0; i<nge; i++)
  {
    xge[i]  = 6 + i;
    yge[i]  = 0;
    eyge[i]  = 3.853e-4;
  }

  TGraphErrors *gr_lum = new TGraphErrors(nge, xge, yge, 0, eyge);
  gr_lum->SetFillColor(4);
  gr_lum->SetLineWidth(1504);
  gr_lum->SetFillStyle(3005);

  mc();
  mcd();
  gPad->SetGridy();
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.04);
  cout.precision(4);
  for(int beam=0; beam<3; beam++)
  {
    for(int ipt=7; ipt<npT_pol-1; ipt++)
    {
      double xpt, comb[3], ecomb[3];
      for(int isys=0; isys<3; isys++)
      {
        int igr = beam + ngr_photon*2 + ngr_photon*3*isys;
        qt_all->Query(ipt, igr, xpt, comb[isys], ecomb[isys]);
      }
      // comb[3] = {"A_{LL}^{default}", "A_{LL}^{#eta}", "#pi^{0} fitting"};
      double esys = sqrt(pow(comb[1]-comb[0],2) + pow(comb[2]-comb[0],2));
      qt_sys->Fill(ipt, beam, xpt, comb[0], esys);
      if( beam==2 && TMath::Finite(comb[0]+ecomb[0]+esys) )
        cout << fixed << xpt << " & " << scientific << comb[0] << " & " << ecomb[0] << " (" << fixed << 100.*ecomb[0]/fabs(comb[0]) << "\\%) & "
          << scientific << esys << " (" << fixed << 100.*esys/fabs(ecomb[0]) << "\\%) \\\\" << endl;
    } // ipt

    int igr = beam + ngr_photon*2;
    TGraphErrors *gr_all = qt_all->Graph(igr);
    TGraphErrors *gr_sys = qt_sys->Graph(beam);
    while(true)
    {
      double xgr, ygr;
      gr_all->GetPoint(0, xgr, ygr);
      if(xgr > 6.) break;
      gr_all->RemovePoint(0);
    }
    gr_all->SetTitle("");
    aset(gr_all, "p_{T} [GeV/c]",beam_list[beam], 4.9,20.1, -0.06,0.05);
    style(gr_all, 20, 1);
    style(gr_sys, 1, 1);
    gr_all->SetMarkerSize(0.8);
    gr_all->GetXaxis()->SetNdivisions(505);
    gr_all->GetYaxis()->SetNdivisions(510);
    gr_sys->SetLineWidth(2);
    gr_all->Draw("AP");
    //gr_lum->Draw("3");
    //gr_sys->Draw("[]");
    for(int i=0; i<gr_sys->GetN(); i++)
    {
      double xx, yy;
      gr_sys->GetPoint(i, xx, yy);
      double eyy = gr_sys->GetErrorY(i);
      box->DrawBox(xx-0.2,yy-eyy,xx+0.2,yy+eyy);
    }
    if(beam==2)
    {
      gr_dssv->SetLineColor(kRed);
      gr_dssv->Draw("C");
      latex->DrawLatexNDC(0.23,0.82, "#splitline{Isolated direct photon A_{LL}}{#vec{p}+#vec{p} #sqrt{s} = 510 GeV, |#eta| < 0.25}");
      latex->DrawLatexNDC(0.23,0.47, "#scale[0.8]{#splitline{3.9e-4 shift uncertainty from}{relative luminosity not included}}");
      latex->DrawLatexNDC(0.23,0.38, "#scale[0.8]{#splitline{6.6% scale uncertainty from}{polarization not included}}");
      latex->DrawLatexNDC(0.75,0.65, "DSSV14");
    }
    const char *outfile = Form("plots/IsoPhotonALL-beam%d", beam);
    c0->Print(Form("%s.pdf", outfile));
    c0->Clear("D");
    if(beam == 2)
    {
      char *cmd = Form("preliminary.pl --input=%s.pdf --output=%s-prelim.pdf --x=150 --y=130 --scale=0.8", outfile,outfile);
      system(cmd);
    }
  } // beam

  qt_sys->Save();
}
