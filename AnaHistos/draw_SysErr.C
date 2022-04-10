#include "GlobalVars.h"
#include "QueryTree.h"
#include "Chi2Fit.h"

void draw_SysErr(const int pwhg = 0, const int ipwhg = 0, const int prelim = 0)
{
  const double PI = TMath::Pi();
  const double DeltaEta = 0.5;

  if(pwhg == 0)
  {
    const char *suffix = "-jetphox";
    const double jetphox_scale = 1./2000.;  // combined 2000 histograms
    const int nmu[2] = {3, 5};
    const char *jetphox_fname[2][nmu[1]] = {
      {"onept", "halfpt", "twopt"},
      {"MMM", "LLH", "LHH", "HLL", "HHL"}
    };
    const char *mu_name[2][3] = {
      {"#mu = p_{T}", "#mu = p_{T}/2", "#mu = 2p_{T}"},
      {"#mu_{R} = 0.56p_{T}, #mu_{F} = p_{T}", "#mu_{R} = 0.54p_{T}, #mu_{F} = p_{T}/2 or 2p_{T}", "#mu_{R} = 0.58p_{T}, #mu_{F} = p_{T}/2 or 2p_{T}"}
    };
  }
  if(pwhg == 1 || pwhg == 3)
  {
    const char *pwhg_type[4] = {"with MPI", "QED-QCD veto", "without MPI", "pure hard"};
    const char *pwhg_suffix[4] = {"-pwhg", "-qedqcd", "-nompi", "-purehard"};
    TFile *f_pythia = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaPowheg-histo%s.root",ipwhg?pwhg_suffix[ipwhg]:""));
    TH1 *h_events = (TH1*)f_pythia->Get("h_events");
    const double nEvents = h_events->GetBinContent(1);
    const int nmu[2] = {7, 7};
    const char *mu_name[3] = {"#mu_{R} = p_{T}, #mu_{F} = p_{T}", "#mu_{R} = p_{T}/2, #mu_{F} = p_{T}/2 or p_{T}", "#mu_{R} = 2p_{T}, #mu_{F} = p_{T} or 2p_{T}"};
  }
  if(pwhg == 2 || pwhg == 3)
  {
    const char *suffix = "-werner";
    const int nmu[2] = {3, 3};
    const char *scale_name[3] = {"nnpdf-grv-onept", "nnpdf-grv-halfpt", "nnpdf-grv-twopt"};
    const char *pwhg_mu_name[3] = {"#mu = p_{T}", "#mu = p_{T}/2", "#mu = 2p_{T}"};
    TGraph *gr_werner[3];
  }
  if(pwhg == 3)
  {
    const char *suffix = "-werner-pwhg";
    TGraphAsymmErrors *gr_band = new TGraphAsymmErrors(npT);
    TGraphAsymmErrors *gr_band_ratio = new TGraphAsymmErrors(npT);
    TGraph *gr_central_ratio = new TGraph(npT);
  }
  if(pwhg < 0 || pwhg > 3)
  {
    cout << "Wrong input" << endl;
    return;
  }

  QueryTree *qt_sys = new QueryTree("data/CrossSection-syserr.root", "RECREATE");
  TGraphErrors *gr_cross[2], *gr_cross_sys[2];

  for(int iso=0; iso<2; iso++)
  {
    QueryTree *qt_cross = new QueryTree(Form("data/CrossSection-%sphoton.root",iso?"iso":""));

    cout.precision(4);
    for(int ipt=12; ipt<npT; ipt++)
    {
      double xpt, xsec, exsec, rsys, ersys;
      qt_cross->Query(ipt, 3, xpt, xsec, exsec);
      qt_cross->Query(ipt, 4, xpt, rsys, ersys);
      double sys = xsec*rsys;
      qt_sys->Fill(ipt, iso, xpt, xsec, sys);
      if( pwhg == 0 && TMath::Finite(xsec+exsec+sys) && xsec > 0. )
      {
        //cout << fixed << xpt << " & " << xsec << " & " << exsec << " (" << setfill('0') << setw(7) << 100.*exsec/xsec << "\\%) & "
        //  << sys << " (" << setfill('0') << setw(7) << 100.*sys/xsec << "\\%) \\\\" << endl;
        cout << fixed << setprecision(1) << pTbin[ipt] << "\t" << pTbin[ipt+1] << "\t"
          << scientific << setprecision(3) << xsec << "\t" << setprecision(1) << exsec << "\t" << sys << endl;
      }
    }
    cout << "***" << endl;

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
    legi(0, 0.25,0.03,0.47,0.20);
    leg0->SetTextSize(0.030);
    legi(1, 0.60,0.35,0.80,0.55);
    leg1->SetTextSize(0.045);
    if(pwhg == 3)
    {
      legi(2, 0.65,0.30,0.90,0.35);
      leg2->SetTextSize(0.030);
    }
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);

    pad2->cd();
    gPad->SetLeftMargin(0.2);
    gPad->SetTopMargin(0.);
    gPad->SetBottomMargin(0.25);

    TLine *line = new TLine();
    line->SetLineWidth(2);

    int igr_band = 0;
    double nlo[npT][nmu[1]], band[npT][nmu[1]];
    for(int imu=0; imu<nmu[iso]; imu++)
    {
      if(pwhg == 0)
      {
        TFile *f_nlo = new TFile( Form("data/%sprompt-x2000-ct14-%s.root",iso?"iso":"inc",jetphox_fname[iso][imu]) );
        TH1 *h_nlo = (TH1*)f_nlo->Get("hp41");
        h_nlo->Scale(jetphox_scale);
      }
      if(pwhg == 1 || pwhg == 3)
      {
        TH1 *h_nlo = (TH1*)f_pythia->Get(Form("hard0_iso%d_rap0_id%d",iso,imu));
        h_nlo->Scale(1./nEvents, "width");
      }
      if(pwhg == 2 || pwhg == 3)
      {
        gr_werner[imu] = new TGraph(Form("data/werner-cross-%s-%s.txt",iso?"iso":"inc",scale_name[imu]));
      }
      for(int ipt=12; ipt<npT; ipt++)
      {
        if(pwhg != 2)
        {
          double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
          int bin_th = h_nlo->GetXaxis()->FindBin(xpt);
          nlo[ipt][imu] = h_nlo->GetBinContent(bin_th);
        }
        if(pwhg == 2 || pwhg == 3)
        {
          double xpt;
          gr_werner[imu]->GetPoint(ipt-4, xpt, band[ipt][imu]);
          if(pwhg == 3 && imu == 2)
          {
            gr_band->SetPoint(igr_band, xpt, band[ipt][0]);
            gr_band->SetPointError(igr_band, 0, 0, band[ipt][0]-band[ipt][2], band[ipt][1]-band[ipt][0]);
            gr_band_ratio->SetPoint(igr_band, xpt, 0);
            gr_band_ratio->SetPointError(igr_band, 0, 0, 1-band[ipt][2]/band[ipt][0], band[ipt][1]/band[ipt][0]-1);
            igr_band++;
          }
        }
      }
      if(pwhg == 0)
        delete f_nlo;
    } // imu
    if(pwhg == 3)
    {
      gr_band->Set(igr_band);
      gr_band_ratio->Set(igr_band);
    }

    TGraphErrors *gr_central;
    for(int imu=0; imu<3; imu++)
    {
      TGraph *gr_nlo = new TGraph(npT);
      TGraphErrors *gr_ratio = new TGraphErrors(npT);
      TGraphErrors *gr_ratio_sys = new TGraphErrors(npT);

      int igp = 0;
      for(int ipt=12; ipt<npT; ipt++)
      {
        double xpt, Combine, eCombine, sysCombine;
        if( !qt_cross->Query(ipt, 3, xpt, Combine, eCombine) ||
            !qt_sys->Query(ipt, iso, xpt, Combine, sysCombine) )
          continue;

        if(pwhg != 2 && imu == 0)
          for(int i=0; i<nmu[iso]; i++)
            nlo[ipt][i] /= (2*PI*xpt*DeltaEta);

        double sigma_nlo;
        if(imu == 0)
          sigma_nlo = pwhg!=2?nlo[ipt][imu]:band[ipt][imu];
        else if(imu == 1)
          sigma_nlo = TMath::MaxElement(nmu[iso], pwhg!=2?nlo[ipt]:band[ipt]);
        else if(imu == 2)
          sigma_nlo = TMath::MinElement(nmu[iso], pwhg!=2?nlo[ipt]:band[ipt]);
        gr_nlo->SetPoint(igp, xpt, sigma_nlo);

        double den = pwhg<2?nlo[ipt][0]:band[ipt][0];
        double ratio, eratio, sysratio;
        if(imu == 0)
        {
          ratio = Combine/den;
          eratio = eCombine/den;
          sysratio = sysCombine/den;
          gr_ratio_sys->SetPoint(igp, xpt, ratio-1);
          gr_ratio_sys->SetPointError(igp, 0., sysratio);
          if(pwhg == 3)
            gr_central_ratio->SetPoint(igp, xpt, nlo[ipt][0]/band[ipt][0]-1);
        }
        else
        {
          ratio = sigma_nlo/den;
          eratio = 1e-9;
        }
        gr_ratio->SetPoint(igp, xpt, ratio-1);
        gr_ratio->SetPointError(igp, 0., eratio);
        igp++;
      } // ipt

      gr_nlo->Set(igp);
      gr_ratio->Set(igp);
      gr_ratio_sys->Set(igp);
      if(pwhg == 3 && imu == 0)
      {
        gr_central_ratio->Set(igp);
        style(gr_central_ratio, 1, 1, 2);
      }

      style(gr_nlo, imu+1, 1, 2);
      style(gr_ratio, imu==0?20:imu+1, imu==0?2:1, 2);
      gr_ratio->SetMarkerSize(0.8);
      if(pwhg != 2)
      {
        if(imu == 0)
          gr_central == gr_nlo;
        else
          leg0->AddEntry(gr_nlo, pwhg==0?mu_name[iso][imu]:mu_name[imu], "L");
        if(imu == 1)
          leg0->AddEntry(gr_central, pwhg==0?mu_name[iso][0]:mu_name[0], "L");
      }

      if(imu == 0)
      {
        pad1->cd();
        if(pwhg == 3)
        {
          aset(gr_band, "p_{T} (GeV/c)", "Ed^{3}#sigma/dp^{3} (pb GeV^{-2} c^{3})", 4.9,30.1, 0.5e-1, iso?2e3:5e3);
          style(gr_band, 1, 1);
          gr_band->SetTitle("");
          gr_band->SetFillColor(kCyan-7);
          //gr_band->SetFillStyle(3001);
          gr_band->Draw("A3");
          gr_band->Draw("X");
          leg2->AddEntry(gr_band, "#mu = p_{T}/2, p_{T}, 2p_{T}", "F");
          leg2->Draw();
        }
        gr_cross[iso] = qt_cross->Graph(3);
        gr_cross_sys[iso] = qt_sys->Graph(iso);
        gr_cross[iso]->SetTitle("");
        aset(gr_cross[iso], "p_{T} (GeV/c)", "Ed^{3}#sigma/dp^{3} (pb GeV^{-2} c^{3})", 4.9,30.1, 0.5e-1, iso?2e3:5e3);
        style(gr_cross[iso], 20, 2, 2);
        style(gr_cross_sys[iso], 1, 2, 2);
        gr_cross[iso]->SetMarkerSize(0.8);
        gr_cross[iso]->Draw(pwhg==3?"P":"AP");
        gr_cross_sys[iso]->Draw("[]");
        gr_nlo->Draw("LX");
        leg0->Draw();
        latex->DrawLatexNDC(0.29,0.87, Form("#splitline{%s direct photon cross section}{p+p #sqrt{s} = 510 GeV, |#eta| < 0.25}",iso?"Isolated":"Inclusive"));
        if(pwhg != 2)
          latex->DrawLatexNDC(0.29,0.79, "#scale[0.8]{10% absolute luminosity uncertainty not shown}");
        if(prelim == 0 && pwhg != 3)
          leg1->AddEntry(gr_cross[iso], "PHENIX Data", "P");
        leg1->Draw();
        if(pwhg == 0)
        {
          latex->DrawLatexNDC(0.25,0.37, "NLO pQCD");
          latex->DrawLatexNDC(0.25,0.32, "(by JETPHOX)");
          latex->DrawLatexNDC(0.25,0.27, "CT14 PDF");
          latex->DrawLatexNDC(0.25,0.22, "BFG II FF");
        }
        if(pwhg == 1 || pwhg == 3)
        {
          latex->DrawLatexNDC(0.25,0.37, "NLO pQCD");
          latex->DrawLatexNDC(0.25,0.32, "(by POWHEG");
          latex->DrawLatexNDC(0.25,0.27, Form("%s)",pwhg_type[ipwhg]));
          latex->DrawLatexNDC(0.25,0.22, "CT14 PDF");
        }
        if(pwhg == 2)
        {
          latex->DrawLatexNDC(0.86,0.91, Form("(%s)",iso?"c":"a"));
          latex->DrawLatexNDC(0.25,0.27, "NLO pQCD");
          latex->DrawLatexNDC(0.25,0.22, "(by W. Vogelsang)");
          latex->DrawLatexNDC(0.25,0.17, "NNPDF3.0 PDF");
          latex->DrawLatexNDC(0.25,0.12, "GRV FF");
          latex->DrawLatexNDC(0.25,0.07, "#mu = p_{T}/2, p_{T}, 2p_{T}");
        }
        if(pwhg == 3)
        {
          latex->DrawLatexNDC(0.60,0.50, "NLO pQCD");
          latex->DrawLatexNDC(0.60,0.45, "(by W. Vogelsang)");
          latex->DrawLatexNDC(0.60,0.40, "NNPDF3.0 PDF");
          latex->DrawLatexNDC(0.60,0.35, "GRV FF");
          if(prelim == 0)
            latex->DrawLatexNDC(0.23,0.50, "#splitline{PHENIX}{Data}");
        }
        if(iso)
        {
          latex->DrawLatexNDC(0.45,0.70, "Isolation cut condition");
          latex->DrawLatexNDC(0.45,0.60, "#splitline{r_{cone} = #sqrt{(#delta#eta)^{2} + (#delta#phi)^{2}} = 0.5}{E_{cone} < 0.1E_{#gamma}}");
        }

        pad2->cd();
        if(pwhg == 3)
        {
          gr_band_ratio->SetTitle(";p_{T} (GeV/c);#scale[0.9]{#frac{Data-Theory}{Theory}}");
          aset(gr_band_ratio, "","", 4.9,30.1, iso?-0.25:-0.45,(iso&&pwhg!=1)?0.45+0.125*pwhg:2.15-iso, 1.,0.6,0.1,0.12);
          gr_band_ratio->GetXaxis()->SetLabelSize(0.09);
          gr_band_ratio->GetYaxis()->SetLabelSize(0.09);
          gr_band_ratio->GetYaxis()->SetLabelOffset(0.005);
          gr_band_ratio->GetXaxis()->SetTickSize(0.08);
          gr_band_ratio->SetFillColor(kCyan-7);
          //gr_band_ratio->SetFillStyle(3001);
          gr_band_ratio->Draw("A3");
          gr_band_ratio->Draw("X");
          gr_central_ratio->Draw("LX");
        }
        gr_ratio->SetTitle(";p_{T} (GeV/c);#scale[0.9]{#frac{Data-Theory}{Theory}}");
        aset(gr_ratio, "","", 4.9,30.1, iso?-0.25:-0.45,(iso&&pwhg!=1)?0.45+0.125*pwhg:2.15-iso, 1.,0.6,0.1,0.12);
        style(gr_ratio_sys, 1, 2, 2);
        gr_ratio->GetXaxis()->SetLabelSize(0.09);
        gr_ratio->GetYaxis()->SetLabelSize(0.09);
        gr_ratio->GetYaxis()->SetLabelOffset(0.005);
        gr_ratio->GetXaxis()->SetTickSize(0.08);
        gr_ratio->Draw(pwhg==3?"PE":"APE");
        gr_ratio_sys->Draw("[]");
        line->DrawLine(pwhg==3?5.:6., 0., 30., 0.);
      }
      else
      {
        pad1->cd();
        gr_nlo->Draw("LX");

        pad2->cd();
        gr_ratio->Draw("LX");
      }
      if(pwhg == 2)
        latex->DrawLatexNDC(0.86,0.93, Form("#scale[1.8]{(%s)}",iso?"d":"b"));
    } // imu

    const char *outfile = Form("plots/CrossSection-%sphoton%s", iso?"iso":"",pwhg==1?pwhg_suffix[ipwhg]:suffix);
    c0->Print(Form("%s.pdf", outfile));
    if(prelim == 1)
    {
      char *cmd = Form("preliminary.pl --input=%s.pdf --output=%s-prelim.pdf --x=360 --y=420 --scale=0.8", outfile,outfile);
      system(cmd);
      cmd = Form("rm %s.pdf", outfile);
      system(cmd);
    }
    delete c0;
  } // iso

  qt_sys->Save();

  if(false)
  {
    TCanvas *c1 = new TCanvas("c1", "c1", 600,800);
    TPad *pad3 = new TPad("pad3", "pad3", 0.,0.,1.,1., -1,0);
    pad3->Draw();

    pad3->cd();
    gPad->SetLeftMargin(0.2);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetLogy();
    legi(2, 0.22,0.18,0.45,0.28);
    leg0->SetTextSize(0.035);
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);

    for(int iso=0; iso<2; iso++)
    {
      style(gr_cross[iso], 20+iso, 1+iso, 2);
      style(gr_cross_sys[iso], 1, 1+iso, 2);
      gr_cross[iso]->SetMarkerSize(0.8);
      gr_cross[iso]->Draw(iso?"P":"AP");
      gr_cross_sys[iso]->Draw("[]");
      leg2->AddEntry(gr_cross[iso], Form("%s direct photon",iso?"Isolated":"Inclusive"), "P");
    }
    leg2->Draw();
    latex->DrawLatexNDC(0.35,0.88, "#splitline{Direct photon cross section}{p+p #sqrt{s} = 510 GeV, |#eta| < 0.25}");
    latex->DrawLatexNDC(0.35,0.80, "#scale[0.8]{#splitline{10% absolute luminosity}{uncertainty not shown}}");
    latex->DrawLatexNDC(0.45,0.70, "Isolation cut condition");
    latex->DrawLatexNDC(0.45,0.62, "#splitline{r_{cone} = #sqrt{(#delta#eta)^{2} + (#delta#phi)^{2}} = 0.5}{E_{cone} < 0.1E_{#gamma}}");

    const char *outfile = "plots/CrossSection-photon-isophoton";
    c1->Print(Form("%s.pdf", outfile));
    char *cmd = Form("preliminary.pl --input=%s.pdf --output=%s-prelim.pdf --x=360 --y=320 --scale=0.8", outfile,outfile);
    system(cmd);
    cmd = Form("rm %s.pdf", outfile);
    system(cmd);
  }
}
