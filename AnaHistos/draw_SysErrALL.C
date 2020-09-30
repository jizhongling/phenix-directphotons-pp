#include "GlobalVars.h"
#include "QueryTree.h"
#include "IsoPhotonALL.h"

void draw_SysErrALL()
{
  const char *beam_list[3] = {"A_{L}^{Blue}", "A_{L}^{Yellow}", "A_{LL}"};
  const char *err_name[4] = {"Rel Lum", "A_{LL}^{#eta}", "#pi^{0} fitting", "Pol"};

  QueryTree *qt_sys = new QueryTree("data/IsoPhotonALL-syserr.root", "RECREATE");
  QueryTree *qt_all = new QueryTree("data/IsoPhotonALL.root");

  TGraph *gr_dssv = new TGraph("data/dssv-all.txt", "%lg %lg");

  TBox *box = new TBox();
  box->SetLineColor(2);
  box->SetFillStyle(0);

  TGraph *gr_ratio[3];
  for(int i=0; i<3; i++)
    gr_ratio[i] = new TGraph(npT_pol);

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
        if(beam==2 && isys>0)
          gr_ratio[isys]->SetPoint(ipt, xpt, fabs(comb[isys]/comb[0]-1));
      }
      double esys = sqrt(pow(3.853e-4,2) + pow(comb[0]*0.066/(beam<2?2:1),2) + pow(comb[2]-comb[0],2) + pow(comb[1]-comb[0],2));
      qt_sys->Fill(ipt, beam, xpt, comb[0], esys);
      if(beam==2)
        gr_ratio[0]->SetPoint(ipt, xpt, 3.853e-4/comb[0]);
      if( beam==2 && TMath::Finite(comb[0]+ecomb[0]+esys) )
        cout << fixed << xpt << " & " << scientific << comb[0] << " & " << ecomb[0] << " (" << fixed << 100.*ecomb[0]/fabs(comb[0]) << "\\%) & "
          << scientific << esys << " (" << fixed << 100.*esys/fabs(ecomb[0]) << "\\%) \\\\" << endl;
    } // ipt

    int igr = beam + ngr_photon*2;
    TGraphErrors *gr_all = qt_all->Graph(igr);
    TGraphErrors *gr_sys = qt_sys->Graph(beam);
    gr_all->SetTitle( Form("#gamma^{dir} %s",beam_list[beam]) );
    aset(gr_all, "p_{T} [GeV/c]",beam_list[beam], 5.9,20.1, -0.06,0.05);
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
    if(beam==2)
    {
      gr_dssv->SetLineColor(kRed);
      gr_dssv->Draw("C");
      latex->DrawLatexNDC(0.23,0.37, "#splitline{Isolated direct photon A_{LL}}{#vec{p}+#vec{p} #sqrt{s} = 510 GeV, |#eta| < 0.25}");
      latex->DrawLatexNDC(0.23,0.27, "#scale[0.8]{#DeltaA_{LL}^{RelLum} = 3.853e-4, #frac{#DeltaA_{LL}^{Pol}}{|A_{LL}|} = 6.6%}");
      latex->DrawLatexNDC(0.75,0.65, "DSSV14");
    }
    c0->Print(Form("plots/IsoPhotonALL-beam%d.pdf",beam));
    c0->Clear("D");
  } // beam
  system("preliminary.pl --input=plots/IsoPhotonALL-beam2.pdf --output=plots/IsoPhotonALL-beam2-prelim.pdf --x=360 --y=420 --scale=0.8");

  mc(1);
  mcd(1);
  legi(0, 0.4,0.7,0.9,0.9);
  leg0->SetNColumns(2);
  for(int i=0; i<3; i++)
  {
    aset(gr_ratio[i], "p_{T} [GeV/c]","Rel Err", 5.9,20.1, 0.,0.15);
    style(gr_ratio[i], 1, 1+i);
    gr_ratio[i]->Draw(i==0?"AL":"L");
    leg0->AddEntry(gr_ratio[i], err_name[i], "L");
  }
  TLine *l_ratio = new TLine;
  l_ratio->SetLineColor(4);
  l_ratio->DrawLine(6.,0.066,17.,0.066);
  leg0->AddEntry(l_ratio, err_name[3], "L");
  leg0->Draw();
  //c1->Print("plots/IsoPhotonALL-beam2-relerr.pdf");

  qt_sys->Save();
}
