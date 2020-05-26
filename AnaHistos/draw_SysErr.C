#include "GlobalVars.h"
#include "QueryTree.h"
#include "Chi2Fit.h"

void draw_SysErr()
{
  const double PI = TMath::Pi();
  const double DeltaEta = 0.5;
  const double jetphox_scale = 1./400.;  // combined 400 histograms
  const char *jetphox_fname[3] = {"onept", "halfpt", "twopt"};
  const char *mu_name[3] = {"#mu=p_{T}", "#mu=p_{T}/2", "#mu=2p_{T}"};

  QueryTree *qt_sys = new QueryTree("data/CrossSection-syserr.root", "RECREATE");

  for(int iso=0; iso<2; iso++)
  {
    char *type = iso ? "iso" : "";
    QueryTree *qt_cross = new QueryTree(Form("data/CrossSection-%sphoton.root",type));

    cout.precision(4);
    for(int ipt=0; ipt<npT; ipt++)
    {
      double xpt, xsec, exsec, rsys, ersys;
      qt_cross->Query(ipt, 3, xpt, xsec, exsec);
      qt_cross->Query(ipt, 4, xpt, rsys, ersys);
      double sys = xsec*rsys;
      qt_sys->Fill(ipt, iso, xpt, xsec, sys);
      if( TMath::Finite(xsec+exsec+sys) && xsec > 0. )
        cout << xpt << " & " << xsec << " & " << exsec << " (" << 100.*exsec/xsec << "\\%) & "
          << sys << " (" << 100.*sys/xsec << "\\%) \\\\" << endl;
    }

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
    legi(0, 0.25,0.05,0.5,0.3);
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);

    pad2->cd();
    gPad->SetLeftMargin(0.2);
    gPad->SetTopMargin(0.);
    gPad->SetBottomMargin(0.25);

    TLine *line = new TLine();
    line->SetLineWidth(2);

    double sigma_onept[npT], esigma_onept[npT];
    for(int imu=0; imu<3; imu++)
    {
      TGraphErrors *gr_nlo = new TGraphErrors(npT);
      TGraphErrors *gr_ratio = new TGraphErrors(npT);
      TGraphErrors *gr_ratio_sys = new TGraphErrors(npT);
      char *type = iso ? "iso" : "inc";
      TFile *f_nlo = new TFile( Form("data/%sprompt-x400-ct14-%s.root",type,jetphox_fname[imu]) );
      TH1 *h_nlo = (TH1*)f_nlo->Get("hp41");
      h_nlo->Scale(jetphox_scale);

      int igp = 0;
      for(int ipt=12; ipt<npT; ipt++)
      {
        double xpt, Combine, eCombine, sysCombine;
        if( !qt_cross->Query(ipt, 3, xpt, Combine, eCombine) ||
            !qt_sys->Query(ipt, iso, xpt, Combine, sysCombine) )
          continue;

        double factor = 1. / (2*PI*xpt*DeltaEta);
        int bin_th = h_nlo->GetXaxis()->FindBin(xpt);
        double sigma_nlo = factor * h_nlo->GetBinContent(bin_th);
        double esigma_nlo = factor * h_nlo->GetBinError(bin_th);
        gr_nlo->SetPoint(igp, xpt, sigma_nlo);
        gr_nlo->SetPointError(igp, 0., esigma_nlo);

        double ratio, eratio, sysratio;
        if(imu == 0)
        {
          sigma_onept[ipt] = sigma_nlo;
          esigma_onept[ipt] = esigma_nlo;
          ratio = Combine/sigma_nlo;
          eratio = ratio*sqrt(pow(eCombine/Combine,2) + pow(esigma_nlo/sigma_nlo,2));
          sysratio = sysCombine/sigma_nlo; 
          gr_ratio_sys->SetPoint(igp, xpt, ratio-1);
          gr_ratio_sys->SetPointError(igp, 0., sysratio);
        }
        else
        {
          ratio = sigma_nlo/sigma_onept[ipt];
          eratio = ratio*sqrt(pow(esigma_onept[ipt]/sigma_onept[ipt],2) + pow(esigma_nlo/sigma_nlo,2));
        }
        gr_ratio->SetPoint(igp, xpt, ratio-1);
        gr_ratio->SetPointError(igp, 0., eratio);
        igp++;
      } // ipt

      gr_nlo->Set(igp);
      gr_ratio->Set(igp);
      gr_ratio_sys->Set(igp);

      style(gr_nlo, 1, imu+1, 2);
      style(gr_ratio, 20, imu+1, 2);
      gr_ratio->SetMarkerSize(0.8);
      leg0->AddEntry(gr_nlo, Form("%s",mu_name[imu]), "L");

      if(imu == 0)
      {
        pad1->cd();
        TGraphErrors *gr_cross = qt_cross->Graph(3);
        TGraphErrors *gr_cross_sys = qt_sys->Graph(iso);
        gr_cross->SetTitle("");
        aset(gr_cross, "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{3}]", 5.9,30.1, 0.5e-1, 2e3);
        style(gr_cross, 20, 1, 2);
        style(gr_cross_sys, 1, 1, 2);
        gr_cross->SetMarkerSize(0.8);
        gr_cross->Draw("AP");
        gr_cross_sys->Draw("[]");
        gr_nlo->Draw("LX");
        leg0->Draw();
        latex->DrawLatexNDC(0.29,0.87, Form("#splitline{%s direct photon cross section}{p+p #sqrt{s} = 510 GeV, |#eta| < 0.25}",iso?"Isolated":"Inclusive"));
        latex->DrawLatexNDC(0.29,0.79, "#scale[0.8]{10% absolute luminosity uncertainty not included}");
        latex->DrawLatexNDC(0.25,0.36, "#splitline{NLO pQCD}{(by JETPHOX)}");
        latex->DrawLatexNDC(0.25,0.30, "CT14 PDF");
        if(iso)
        {
          latex->DrawLatexNDC(0.45,0.70, "Isolation cut condition");
          latex->DrawLatexNDC(0.45,0.60, "#splitline{#sqrt{(#Delta#eta)^{2} + (#Delta#phi)^{2}} < 0.5}{E_{cone} < 0.1E_{#gamma}}");
        }

        pad2->cd();
        gr_ratio->SetTitle(";p_{T} [GeV/c];#frac{Data-Theory}{Theory}");
        aset(gr_ratio, "","", 5.9,30.1, iso?-0.23:-0.63,iso?0.63:3.13, 1.,0.6,0.1,0.12);
        style(gr_ratio_sys, 1, imu+1, 2);
        gr_ratio->GetXaxis()->SetLabelSize(0.08);
        gr_ratio->GetYaxis()->SetLabelSize(0.08);
        gr_ratio->GetXaxis()->SetTickSize(0.08);
        gr_ratio->Draw("APE");
        gr_ratio_sys->Draw("[]");
        line->DrawLine(6., 0., 30., 0.);
      }
      else
      {
        pad1->cd();
        gr_nlo->Draw("LX");

        pad2->cd();
        gr_ratio->Draw("LX");
      }

      delete h_nlo;
      delete f_nlo;
    }

    char *type = iso ? "iso" : "";
    c0->Print(Form("plots/CrossSection-%sphoton-syserr.pdf",type));
    const char *cmd = Form("preliminary.pl --input=plots/CrossSection-%sphoton-syserr.pdf --output=plots/CrossSection-%sphoton-prelim.pdf --x=360 --y=420 --scale=0.8", type,type);
    system(cmd);
    delete c0;
  }

  qt_sys->Save();
}
