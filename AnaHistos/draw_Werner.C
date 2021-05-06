#include "GlobalVars.h"
#include "QueryTree.h"
#include "Chi2Fit.h"

void draw_Werner()
{
  const double PI = TMath::Pi();
  const double DeltaEta = 0.5;
  const double jetphox_scale = 1./2000.;  // combined 2000 histograms

  const int nmu[2] = {5, 9};
  const int nana[2] = {4, 8};
  const char *scale_name[2][nmu[1]] = {
    {"nnpdf-grv-onept", "nnpdf-grv-halfpt", "nnpdf-grv-twopt", "ct14-bfg2-onept", "onept"},
    {"nnpdf-grv-onept", "nnpdf-grv-halfpt", "nnpdf-grv-twopt", /*"ct14-grv-onept",*/ "ct14-bfg2-onept", "ct14-bfg2-scR056-scF1", "ct14-grv-scR05-scF05", "ct14-grv-scR05-scF1", "ct14-grv-scR05-scF2", "MMM"}
  };
  const char *leg_name[2][nmu[1]] = {
    {"NNPDF GRV p_{T}", "NNPDF GRV p_{T}/2", "NNPDF GRV 2p_{T}", "CT14 BFGII p_{T}", "JETPHOX p_{T}"},
    {"NNPDF GRV p_{T}", "NNPDF GRV p_{T}/2", "NNPDF GRV 2p_{T}", /*"CT14 GRV p_{T}",*/ "CT14 BFGII p_{T}", "CT14 BFGII #mu_{R}=#mu_{PMC}, #mu_{F}=p_{T}", "CT14 GRV #mu_{R}#approx#mu_{PMC}, #mu_{F}=p_{T}/2", "CT14 GRV #mu_{R}#approx#mu_{PMC}, #mu_{F}=p_{T}", "CT14 GRV #mu_{R}#approx#mu_{PMC}, #mu_{F}=2p_{T}", "JETPHOX #mu_{R}=#mu_{PMC}, #mu_{F}=p_{T}"}
  };

  for(int iso=0; iso<2; iso++)
  {
    TGraph *gr_werner[nmu[1]];
    TGraphErrors *gr_cross, *gr_cross_sys;
    QueryTree *qt_cross = new QueryTree(Form("data/CrossSection-%sphoton.root",iso?"iso":""));
    QueryTree *qt_sys = new QueryTree("data/CrossSection-syserr.root");

    TCanvas *c0 = new TCanvas("c0", "c0", 600,800);
    TPad *pad1 = new TPad("pad1", "pad1", 0.,0.35,1.,1., -1,0);
    TPad *pad2 = new TPad("pad2", "pad2", 0.,0.,1.,0.35, -1,0);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    gPad->SetLeftMargin(0.2);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.);
    gPad->SetLogy();
    legi(0, 0.24-0.02*iso,0.03,0.45,0.33+0.07*iso);
    leg0->SetTextSize(0.035);
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);

    pad2->cd();
    gPad->SetLeftMargin(0.2);
    gPad->SetTopMargin(0.);
    gPad->SetBottomMargin(0.25);

    TLine *line = new TLine();
    line->SetLineWidth(2);

    double nlo[npT][nmu[1]];
    TGraphErrors *gr_central;
    for(int imu=0; imu<nmu[iso]; imu++)
    {
      TFile *f_nlo;
      TH1 *h_nlo;

      if(imu < nana[iso])
      {
        gr_werner[imu] = new TGraph(Form("data/werner-cross-%s-%s.txt",iso?"iso":"inc",scale_name[iso][imu]));
      }
      else
      {
        f_nlo = new TFile( Form("data/%sprompt-x2000-ct14-%s.root",iso?"iso":"inc",scale_name[iso][imu]) );
        h_nlo = (TH1*)f_nlo->Get("hp41");
        h_nlo->Scale(jetphox_scale);
        gr_werner[imu] = new TGraphErrors(npT);
      }

      TGraphErrors *gr_ratio = new TGraphErrors(npT);
      TGraphErrors *gr_ratio_sys = new TGraphErrors(npT);

      int igp = 0;
      for(int ipt=12; ipt<npT; ipt++)
      {
        double xpt, dummy;
        double ratio, eratio, sysratio;
        double Combine, eCombine, sysCombine;

        if( !qt_cross->Query(ipt, 3, xpt, Combine, eCombine) ||
            !qt_sys->Query(ipt, iso, xpt, Combine, sysCombine) )
          continue;

        if(imu < nana[iso])
        {
          gr_werner[imu]->GetPoint(ipt-4, dummy, nlo[ipt][imu]);
        }
        else
        {
          int bin_th = h_nlo->GetXaxis()->FindBin(xpt);
          nlo[ipt][imu] = h_nlo->GetBinContent(bin_th) / (2*PI*xpt*DeltaEta);
          gr_werner[imu]->SetPoint(igp, xpt, nlo[ipt][imu]);
        }

        if(imu == 0)
        {
          ratio = Combine/nlo[ipt][imu];
          eratio = eCombine/nlo[ipt][imu];
          sysratio = sysCombine/nlo[ipt][imu]; 
          gr_ratio_sys->SetPoint(igp, xpt, ratio-1);
          gr_ratio_sys->SetPointError(igp, 0., sysratio);
        }
        else
        {
          ratio = nlo[ipt][imu]/nlo[ipt][0];
          eratio = 1e-9;
        }
        gr_ratio->SetPoint(igp, xpt, ratio-1);
        gr_ratio->SetPointError(igp, 0., eratio);
        igp++;
      } // ipt

      if(imu >= nana[iso])
      {
        gr_werner[imu]->Set(igp);
        delete f_nlo;
      }
      gr_ratio->Set(igp);
      gr_ratio_sys->Set(igp);

      style(gr_werner[imu], imu+1, imu<4||imu>7?imu+1:imu+2, 2);
      style(gr_ratio, imu==0?20:imu+1, imu<4||imu>7?imu+1:imu+2, 2);
      gr_ratio->SetMarkerSize(0.8);
      if(imu == 0)
        gr_central == gr_werner[imu];
      else
        leg0->AddEntry(gr_werner[imu], leg_name[iso][imu], "L");
      if(imu == 1)
        leg0->AddEntry(gr_central, leg_name[iso][0], "L");

      if(imu == 0)
      {
        pad1->cd();
        gr_cross = qt_cross->Graph(3);
        gr_cross_sys = qt_sys->Graph(iso);
        gr_cross->SetTitle("");
        aset(gr_cross, "p_{T} (GeV/c)", "Ed^{3}#sigma/dp^{3} (pb GeV^{-2} c^{3})", 4.9,30.1, 0.5e-1, iso?2e3:5e3);
        style(gr_cross, 20, 1, 2);
        style(gr_cross_sys, 1, 1, 2);
        gr_cross->SetMarkerSize(0.8);
        gr_cross->Draw("AP");
        gr_cross_sys->Draw("[]");
        gr_werner[imu]->Draw("LX");
        latex->DrawLatexNDC(0.29,0.87, Form("#splitline{%s direct photon cross section}{p+p #sqrt{s} = 510 GeV, |#eta| < 0.25}",iso?"Isolated":"Inclusive"));
        latex->DrawLatexNDC(0.29,0.79, "#scale[0.8]{10% absolute luminosity uncertainty not included}");
        latex->DrawLatexNDC(0.25-0.02*iso,0.40+0.07*iso, "NLO pQCD");
        latex->DrawLatexNDC(0.24-0.02*iso,0.35+0.07*iso, "(by W. Vogelsang)");
        if(iso)
        {
          latex->DrawLatexNDC(0.45,0.70, "Isolation cut condition");
          latex->DrawLatexNDC(0.45,0.60, "#splitline{r_{cone} = #sqrt{(#delta#eta)^{2} + (#delta#phi)^{2}} = 0.5}{E_{cone} < 0.1E_{#gamma}}");
        }
        leg0->Draw();

        pad2->cd();
        gr_ratio->SetTitle(";p_{T} (GeV/c);#frac{Data-Theory}{Theory}");
        aset(gr_ratio, "","", 4.9,30.1, iso?-0.25:-0.45,iso?0.70:2.15-iso, 1.,0.6,0.1,0.12);
        style(gr_ratio_sys, 1, imu+1, 2);
        gr_ratio->GetXaxis()->SetLabelSize(0.09);
        gr_ratio->GetYaxis()->SetLabelSize(0.09);
        gr_ratio->GetYaxis()->SetLabelOffset(0.01);
        gr_ratio->GetXaxis()->SetTickSize(0.08);
        gr_ratio->Draw("APE");
        gr_ratio_sys->Draw("[]");
        line->DrawLine(6., 0., 30., 0.);
      }
      else
      {
        pad1->cd();
        gr_werner[imu]->Draw("LX");

        pad2->cd();
        gr_ratio->Draw("LX");
      }
    } // imu

    c0->Print(Form("plots/Werner-%sphoton.pdf",iso?"iso":""));
    delete c0;
  } // iso
}
