#include "GlobalVars.h"
#include "QueryTree.h"
#include "IsoPhotonALL.h"

void draw_SysErrALL(const int prelim = 0)
{
  const char *beam_list[3] = {"A_{L}^{Blue}", "A_{L}^{Yellow}", "A_{LL}"};

  QueryTree *qt_sys = new QueryTree("data/IsoPhotonALL-syserr.root", "RECREATE");
  QueryTree *qt_all = new QueryTree("data/IsoPhotonALL.root");

  TGraphErrors *gr_dssv = new TGraphErrors("data/werner-all-dssv14-nnpdf-grv.txt", "%lg %lg %lg");

  TBox *box = new TBox();
  box->SetLineColor(2);
  box->SetFillStyle(0);

  legi(0, 0.23,0.30,0.50,0.45);
  leg0->SetTextSize(0.030);

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
      {
        //cout << fixed << xpt << " & " << scientific << comb[0] << " & " << ecomb[0] << " (" << fixed << setfill('0') << setw(8) << 100.*ecomb[0]/fabs(comb[0]) << "\\%) & "
        //  << scientific << esys << " (" << fixed << 100.*esys/fabs(ecomb[0]) << "\\%) \\\\" << endl;
        cout << fixed << setprecision(1) << pTbin_pol[ipt] << "\t" << pTbin_pol[ipt+1] << "\t"
          << scientific << setprecision(4) << comb[0] << "\t" << ecomb[0] << "\t" << esys << endl;
      }
    } // ipt
    if(beam == 2)
      cout << "***" << endl;

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
    aset(gr_all, "p_{T} (GeV/c)",beam_list[beam], 4.9,20.1, -0.045,0.05);
    style(gr_all, 20, 2);
    style(gr_sys, 1, 2);
    gr_all->SetMarkerSize(0.8);
    gr_all->GetXaxis()->SetNdivisions(505);
    gr_all->GetYaxis()->SetNdivisions(505);
    gr_sys->SetLineWidth(2);
    gr_all->Draw("AP");
    if(beam==2)
    {
      style(gr_dssv, 1, 1);
      gr_dssv->SetFillColor(kCyan-7);
      //gr_dssv->SetFillStyle(3001);
      gr_dssv->Draw("3");
      gr_dssv->Draw("CX");
      gr_all->Draw("P");
      latex->DrawLatexNDC(0.23,0.85, "#scale[0.9]{#vec{p} + #vec{p} #rightarrow #gamma^{iso} + X, #sqrt{s} = 510 GeV, |#eta| < 0.25}");
      //latex->DrawLatexNDC(0.23,0.79, "#scale[0.6]{3.9#times10^{-4} shift uncertainty from relative luminosity not shown}");
      //latex->DrawLatexNDC(0.23,0.74, "#scale[0.6]{6.6% scale uncertainty from polarization not shown}");
      leg0->AddEntry(gr_all, "PHENIX Data", "P");
      leg0->AddEntry(gr_dssv, "DSSV14 with DSSV_{MC} uncertainty", "LF");
      leg0->Draw();
    }
    //gr_lum->Draw("3");
    //gr_sys->Draw("[]");
    for(int i=0; i<gr_sys->GetN(); i++)
    {
      double xx, yy;
      gr_sys->GetPoint(i, xx, yy);
      double eyy = gr_sys->GetErrorY(i);
      box->DrawBox(xx-0.4,yy-eyy,xx+0.4,yy+eyy);
    }
    const char *outfile = Form("plots/IsoPhotonALL-beam%d", beam);
    c0->Print(Form("%s.pdf", outfile));
    c0->Clear("D");
    if(prelim == 1 && beam == 2)
    {
      char *cmd = Form("preliminary.pl --input=%s.pdf --output=%s-prelim.pdf --x=370 --y=340 --scale=0.8", outfile,outfile);
      system(cmd);
    }
  } // beam

  qt_sys->Save();
}
