#include "GlobalVars.h"
#include "QueryTree.h"
#include "Chi2Fit.h"

void draw_SysErr(const int pwhg = 0, const int ipwhg = 0)
{
  const double PI = TMath::Pi();
  const double DeltaEta = 0.5;
  const char *pwhg_type[4] = {"with MPI", "QED-QCD veto", "without MPI", "pure hard"};
  const char *suffix[4] = {"-pwhg", "-qedqcd", "-nompi", "-purehard"};

  if(pwhg == 0)
  {
    const double jetphox_scale = 1./400.;  // combined 400 histograms
    const int nmu[2] = {3, 7};
    const char *jetphox_fname[2][nmu[1]] = {
      {"onept", "halfpt", "twopt"},
      {"MMM", "LLH", "LHH", "LLL", "HLL", "HHL", "HHH"}
    };
    const char *mu_name[2][3] = {
      {"   1           1            1", " 1/2        1/2         1/2", "   2           2            2"},
      {"0.56         1            1", "0.54    1/2 or 2       2", "0.58    1/2 or 2     1/2"}
    };
  }
  else if(pwhg == 1)
  {
    TFile *f_pythia = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaPowheg-histo%s.root",ipwhg?suffix[ipwhg]:""));
    TH1 *h_events = (TH1*)f_pythia->Get("h_events");
    const double nEvents = h_events->GetBinContent(1);
    const int nmu[2] = {7, 7};
    const char *mu_name[2][3] = {
      {"  1            1           --", "vary       vary        --", "vary       vary        --"},
      {"  1            1           --", "vary       vary        --", "vary       vary        --"}
    };
  }
  else
  {
    cout << "Wrong input" << endl;
    return;
  }

  QueryTree *qt_sys = new QueryTree("data/CrossSection-syserr.root", "RECREATE");
  TGraphErrors *gr_cross[2], *gr_cross_sys[2];

  for(int iso=0; iso<2; iso++)
  {
    char *type = iso ? "iso" : "";
    QueryTree *qt_cross = new QueryTree(Form("data/CrossSection-%sphoton.root",type));

    cout.precision(4);
    for(int ipt=12; ipt<npT; ipt++)
    {
      double xpt, xsec, exsec, rsys, ersys;
      qt_cross->Query(ipt, 3, xpt, xsec, exsec);
      qt_cross->Query(ipt, 4, xpt, rsys, ersys);
      double sys = xsec*rsys;
      qt_sys->Fill(ipt, iso, xpt, xsec, sys);
      if( pwhg == 0 && TMath::Finite(xsec+exsec+sys) && xsec > 0. )
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
    legi(0, 0.22,0.03,0.47,0.20);
    leg0->SetTextSize(0.035);
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);

    pad2->cd();
    gPad->SetLeftMargin(0.2);
    gPad->SetTopMargin(0.);
    gPad->SetBottomMargin(0.25);

    TLine *line = new TLine();
    line->SetLineWidth(2);

    double nlo[npT][nmu[1]], enlo[npT][nmu[1]];
    for(int imu=0; imu<nmu[iso]; imu++)
    {
      char *type = iso ? "iso" : "inc";
      if(pwhg == 0)
      {
        TFile *f_nlo = new TFile( Form("data/%sprompt-x400-ct14-%s.root",type,jetphox_fname[iso][imu]) );
        TH1 *h_nlo = (TH1*)f_nlo->Get("hp41");
        h_nlo->Scale(jetphox_scale);
      }
      else if(pwhg == 1)
      {
        TH1 *h_nlo = (TH1*)f_pythia->Get(Form("hard0_iso%d_rap0_id%d",iso,imu));
        h_nlo->Scale(1./nEvents, "width");
      }
      for(int ipt=12; ipt<npT; ipt++)
      {
        double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
        int bin_th = h_nlo->GetXaxis()->FindBin(xpt);
        nlo[ipt][imu] = h_nlo->GetBinContent(bin_th);
        enlo[ipt][imu] = h_nlo->GetBinError(bin_th);
      }
      if(pwhg == 0)
        delete f_nlo;
    } // imu

    TGraphErrors *gr_central;
    for(int imu=0; imu<3; imu++)
    {
      TGraphErrors *gr_nlo = new TGraphErrors(npT);
      TGraphErrors *gr_ratio = new TGraphErrors(npT);
      TGraphErrors *gr_ratio_sys = new TGraphErrors(npT);

      int igp = 0;
      for(int ipt=12; ipt<npT; ipt++)
      {
        double xpt, Combine, eCombine, sysCombine;
        if( !qt_cross->Query(ipt, 3, xpt, Combine, eCombine) ||
            !qt_sys->Query(ipt, iso, xpt, Combine, sysCombine) )
          continue;

        double factor = 1. / (2*PI*xpt*DeltaEta);
        double sigma_nlo, esigma_nlo;
        if(imu == 0)
        {
          sigma_nlo = nlo[ipt][0];
          esigma_nlo = enlo[ipt][0];
        }
        else if(imu == 1)
        {
          sigma_nlo = TMath::MaxElement(nmu[iso], nlo[ipt]);
          esigma_nlo = TMath::MaxElement(nmu[iso], enlo[ipt]);
        }
        else if(imu == 2)
        {
          sigma_nlo = TMath::MinElement(nmu[iso], nlo[ipt]);
          esigma_nlo = TMath::MaxElement(nmu[iso], enlo[ipt]);
        }
        sigma_nlo *= factor;
        esigma_nlo *= factor;
        gr_nlo->SetPoint(igp, xpt, sigma_nlo);
        gr_nlo->SetPointError(igp, 0., esigma_nlo);

        double ratio, eratio, sysratio;
        if(imu == 0)
        {
          ratio = Combine/sigma_nlo;
          eratio = ratio*sqrt(pow(eCombine/Combine,2) + pow(esigma_nlo/sigma_nlo,2));
          sysratio = sysCombine/sigma_nlo; 
          gr_ratio_sys->SetPoint(igp, xpt, ratio-1);
          gr_ratio_sys->SetPointError(igp, 0., sysratio);
        }
        else
        {
          ratio = sigma_nlo/nlo[ipt][0]/factor;
          eratio = ratio*sqrt(pow(enlo[ipt][0]/nlo[ipt][0],2) + pow(esigma_nlo/sigma_nlo,2));
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
      if(imu == 0)
        gr_central == gr_nlo;
      else
        leg0->AddEntry(gr_nlo, mu_name[iso][imu], "L");
      if(imu == 1)
        leg0->AddEntry(gr_central, mu_name[iso][0], "L");

      if(imu == 0)
      {
        pad1->cd();
        gr_cross[iso] = qt_cross->Graph(3);
        gr_cross_sys[iso] = qt_sys->Graph(iso);
        gr_cross[iso]->SetTitle("");
        aset(gr_cross[iso], "p_{T} [GeV/c]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{3}]", 4.9,30.1, 0.5e-1, iso?2e3:5e3);
        style(gr_cross[iso], 20, 1, 2);
        style(gr_cross_sys[iso], 1, 1, 2);
        gr_cross[iso]->SetMarkerSize(0.8);
        gr_cross[iso]->Draw("AP");
        gr_cross_sys[iso]->Draw("[]");
        gr_nlo->Draw("LX");
        leg0->Draw();
        latex->DrawLatexNDC(0.29,0.87, Form("#splitline{%s direct photon cross section}{p+p #sqrt{s} = 510 GeV, |#eta| < 0.25}",iso?"Isolated":"Inclusive"));
        latex->DrawLatexNDC(0.29,0.79, "#scale[0.8]{10% absolute luminosity uncertainty not included}");
        latex->DrawLatexNDC(0.24,0.40, "NLO pQCD");
        if(pwhg == 0)
        {
          latex->DrawLatexNDC(0.24,0.35, "(by JETPHOX)");
          latex->DrawLatexNDC(0.24,0.28, "#splitline{CT14 PDF}{BFG II FF}");
        }
        else if(pwhg == 1)
        {
          latex->DrawLatexNDC(0.24,0.32, Form("#splitline{(by POWHEG}{%s)}",pwhg_type[ipwhg]));
          latex->DrawLatexNDC(0.24,0.25, "CT14 PDF");
        }
        latex->DrawLatexNDC(0.28,0.22, "#scale[0.8]{#mu_{R}/p_{T}    #mu_{f}/p_{T}    #mu_{F}/p_{T}}");
        if(iso)
        {
          latex->DrawLatexNDC(0.45,0.70, "Isolation cut condition");
          latex->DrawLatexNDC(0.45,0.60, "#splitline{r_{cone} = #sqrt{(#delta#eta)^{2} + (#delta#phi)^{2}} = 0.5}{E_{cone} < 0.1E_{#gamma}}");
        }

        pad2->cd();
        gr_ratio->SetTitle(";p_{T} [GeV/c];#frac{Data-Theory}{Theory}");
        aset(gr_ratio, "","", 4.9,30.1, iso*(1-pwhg)?-0.25:-0.45,iso*(1-pwhg)?0.55:2.15-iso, 1.,0.6,0.1,0.12);
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
        gr_nlo->Draw("LX");

        pad2->cd();
        gr_ratio->Draw("LX");
      }
    } // imu

    char *type = iso ? "iso" : "";
    const char *outfile = Form("plots/CrossSection-%sphoton%s", type,pwhg?suffix[ipwhg]:"-jetphox");
    c0->Print(Form("%s.pdf", outfile));
    if(iso == 0 && ipwhg != 1)
    {
      char *cmd = Form("preliminary.pl --input=%s.pdf --output=%s-prelim.pdf --x=360 --y=420 --scale=0.8", outfile,outfile);
      system(cmd);
      cmd = Form("rm %s.pdf", outfile);
      system(cmd);
    }
    delete c0;
  } // iso

  qt_sys->Save();

  if(pwhg == 0)
  {
    TCanvas *c1 = new TCanvas("c1", "c1", 600,800);
    TPad *pad3 = new TPad("pad3", "pad3", 0.,0.,1.,1., -1,0);
    pad3->Draw();

    pad3->cd();
    gPad->SetLeftMargin(0.2);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetLogy();
    legi(1, 0.22,0.18,0.45,0.28);
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
      leg1->AddEntry(gr_cross[iso], Form("%s direct photon",iso?"Isolated":"Inclusive"), "P");
    }
    leg1->Draw();
    latex->DrawLatexNDC(0.35,0.88, "#splitline{Direct photon cross section}{p+p #sqrt{s} = 510 GeV, |#eta| < 0.25}");
    latex->DrawLatexNDC(0.35,0.80, "#scale[0.8]{#splitline{10% absolute luminosity}{uncertainty not included}}");
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
