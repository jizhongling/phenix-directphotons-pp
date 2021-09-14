#include "GlobalVars.h"
#include "QueryTree.h"
#include "DivideFunctions.h"

void draw_Iso2Inc()
{
  const int sector = 3;  // PbSc west: 0; PbSc east: 1; PbGl: 2; Combined: 3
  const int nmu = 5;
  const int nana = 3;
  const char *scale_name[nmu] = {"nnpdf-grv-halfpt", "nnpdf-grv-onept", "nnpdf-grv-twopt", "", "-nompi"};
  const char *leg_name[nmu] = {"#mu = p_{T}/2", "#mu = p_{T}", "#mu = 2p_{T}", "with MPI", "without MPI"};

  QueryTree *qt_photon = new QueryTree("data/CrossSection-photon.root");
  QueryTree *qt_isophoton = new QueryTree("data/CrossSection-isophoton.root");
  QueryTree *qt_sys = new QueryTree("data/CrossSection-syserr.root");
  QueryTree *qt_pion = new QueryTree("data/CrossSection-pion.root");

  mc();
  mcd();
  legi(0, 0.42,0.25,1.00,0.30);
  legi(1, 0.23,0.67,0.70,0.82);
  legi(2, 0.23,0.20,0.40,0.35);
  leg0->SetTextSize(0.035);
  leg1->SetTextSize(0.035);
  leg2->SetTextSize(0.035);
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.035);

  TGraphAsymmErrors *gr_werner = new TGraphAsymmErrors(npT);
  double ratio_nlo[nana];
  for(int imu=0; imu<nmu; imu++)
  {
    if(imu < nana)
    {
      TGraph *gr_inc = new TGraph(Form("data/werner-cross-inc-%s.txt",scale_name[imu]));
      TGraph *gr_iso = new TGraph(Form("data/werner-cross-iso-%s.txt",scale_name[imu]));

      int igr_nlo = 0;
      for(int ipt=12; ipt<npT; ipt++)
      {
        double xpt, nlo_inc, nlo_iso;
        gr_inc->GetPoint(ipt-4, xpt, nlo_inc);
        gr_iso->GetPoint(ipt-4, xpt, nlo_iso);
        ratio_nlo[imu] = nlo_iso / nlo_inc;
        if( imu==nana-1 && TMath::Finite(ratio_nlo[1]) )
        {
          gr_werner->SetPoint(igr_nlo, xpt, ratio_nlo[1]);
          gr_werner->SetPointError(igr_nlo, 0., 0., ratio_nlo[1]-ratio_nlo[2], ratio_nlo[0]-ratio_nlo[1]);
          igr_nlo++;
        }
      }

      if( imu==nana-1 )
      {
        gr_werner->Set(igr_nlo);
        aset(gr_werner, "p_{T} (GeV/c)", "#gamma_{dir}^{iso}/#gamma_{dir}^{inc}", 4.9,30.1, 0.,1.6);
        style(gr_werner, 1, 1);
        gr_werner->SetTitle("");
        gr_werner->SetFillColor(kCyan-7);
        //gr_werner->SetFillStyle(3001);
        gr_werner->Draw("A3");
        gr_werner->Draw("X");
        leg0->AddEntry(gr_werner, "#mu = p_{T}/2, p_{T}, 2p_{T}", "F");
        latex->DrawLatexNDC(0.43,0.38, "NLO pQCD (by W. Vogelsang)");
        latex->DrawLatexNDC(0.43,0.33, "NNPDF3.0 PDF & GRV FF");

        TGraphErrors *gr[2];
        int igp[2] = {};
        for(int iph=0; iph<2; iph++)
          gr[iph] = new TGraphErrors(npT);
        TGraphErrors *gr_sys = new TGraphErrors(npT);

        for(int ipt=12; ipt<npT; ipt++)
        {
          double xpt, incl, eincl, iso, eiso;

          if( qt_photon->Query(ipt, sector, xpt, incl, eincl) &&
              qt_isophoton->Query(ipt, sector, xpt, iso, eiso) )
          {
            double yy = iso / incl;
            double eyy = yy * sqrt( pow(eincl/incl,2) + pow(eiso/iso,2) );
            if( TMath::Finite(yy+eyy) )
            {
              gr[0]->SetPoint(igp[0], xpt, yy);
              gr[0]->SetPointError(igp[0], 0., eyy);
              double xsec[2], sys[2];
              for(int i=0; i<2; i++)
                qt_sys->Query(ipt, i, xpt, xsec[i], sys[i]);
              double yy = xsec[1] / xsec[0];
              double eyy = yy * sqrt( pow(sys[0]/xsec[0],2) + pow(sys[1]/xsec[1],2) );
              gr_sys->SetPoint(igp[0], xpt, yy);
              gr_sys->SetPointError(igp[0], 0., eyy);
              igp[0]++;
            }
          } // fill photon data

          if( qt_pion->Query(ipt, sector, xpt, incl, eincl) &&
              qt_pion->Query(ipt, sector+4, xpt, iso, eiso) )
          {
            double yy = iso / incl;
            double eyy = yy * sqrt( pow(eincl/incl,2) + pow(eiso/iso,2) );
            if( TMath::Finite(yy+eyy) )
            {
              gr[1]->SetPoint(igp[1], xpt, yy);
              gr[1]->SetPointError(igp[1], 0., eyy);
              igp[1]++;
            }
          } // fill pion data
        } // ipt

        for(int iph=0; iph<1; iph++)
        {
          gr[iph]->Set(igp[iph]);
          style(gr[iph], iph+20, 2);
          if(iph==0)
          {
            gr[iph]->Draw("P");
            style(gr_sys, 1, 2);
            gr_sys->SetLineWidth(2);
            gr_sys->Draw("[]");
            leg2->AddEntry(gr[iph], "#splitline{PHENIX}{Data}", "P");
          } // draw photon ratio
          else
          {
            gr[iph]->Draw("P");
            leg2->AddEntry(gr[iph], "#splitline{PHENIX}{Data #pi^{0}}", "P");
          } // draw pion ratio
        } // iph
      } // imu==nana-1
    } // werner ratio

    else
    {
      TFile *f_pythia = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaPowheg-histo%s.root",scale_name[imu]));
      TH1 *h_photon = (TH1*)f_pythia->Get("hard0_iso0_rap0_id0");
      TH1 *h_isophoton = (TH1*)f_pythia->Get("hard0_iso1_rap0_id0");
      TGraphErrors *gr_pythia = DivideHisto(h_isophoton, h_photon);

      while(true)
      {
        double xgr, ygr;
        gr_pythia->GetPoint(0, xgr, ygr);
        if(xgr > 6.) break;
        gr_pythia->RemovePoint(0);
      }
      for(int i=0; i<gr_pythia->GetN(); i++)
        gr_pythia->SetPointError(i, 0., gr_pythia->GetErrorY(i));

      style(gr_pythia, imu-2, 1);
      gr_pythia->Draw("LE");
      leg1->AddEntry(gr_pythia, leg_name[imu], "L");
      latex->DrawLatexNDC(0.25,0.83, "POWHEG + PYTHIA8");
    } // powheg ratio
  } // imu

  leg0->Draw();
  leg1->Draw();
  leg2->Draw();

  const char *outfile = "plots/Iso2Inc";
  c0->Print(Form("%s.pdf", outfile));
  if(false)
  {
    char *cmd = Form("preliminary.pl --input=%s.pdf --output=%s-prelim.pdf --x=320 --y=200 --scale=0.8", outfile,outfile);
    system(cmd);
    cmd = Form("rm %s.pdf", outfile);
    system(cmd);
  }
}
