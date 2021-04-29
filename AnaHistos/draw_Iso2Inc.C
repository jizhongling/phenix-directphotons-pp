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
  QueryTree *qt_pion = new QueryTree("data/CrossSection-pion.root");
  QueryTree *qt_sys = new QueryTree("data/CrossSection-syserr.root");

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

  mc();
  mcd();
  legi(0, 0.22,0.80,0.90,0.90);
  leg0->SetNColumns(3);
  legi(1, 0.22,0.23,0.90,0.33);
  leg1->SetNColumns(3);
  for(int iph=0; iph<1; iph++)
  {
    gr[iph]->Set(igp[iph]);
    aset(gr[iph], "p_{T} [GeV/c]", "#frac{Isolated}{Inclusive}", 4.9,30.1, 0.,1.6);
    style(gr[iph], iph+20, iph+1);
    if(iph==0)
    {
      gr[iph]->SetTitle("");
      gr[iph]->Draw("AP");
      style(gr_sys, 1, iph+1);
      gr_sys->SetLineWidth(2);
      gr_sys->Draw("[]");
      leg0->AddEntry(gr[iph], "#scale[0.8]{Data}", "P");
    } // draw photon ratio
    else
    {
      gr[iph]->Draw("P");
      leg0->AddEntry(gr[iph], "#scale[0.8]{Data #pi^{0}}", "P");
    } // draw pion ratio
  } // iph

  for(int imu=0; imu<nmu; imu++)
  {
    if(imu < nana)
    {
      TGraph *gr_inc = new TGraph(Form("data/werner-cross-inc-%s.txt",scale_name[imu]));
      TGraph *gr_iso = new TGraph(Form("data/werner-cross-iso-%s.txt",scale_name[imu]));
      TGraph *gr_werner = new TGraphErrors(npT);

      int igr_nlo = 0;
      for(int ipt=12; ipt<npT; ipt++)
      {
        double xpt, nlo_inc, nlo_iso;
        gr_inc->GetPoint(ipt-4, xpt, nlo_inc);
        gr_iso->GetPoint(ipt-4, xpt, nlo_iso);
        double ratio = nlo_iso / nlo_inc;
        if( TMath::Finite(ratio) )
        {
          gr_werner->SetPoint(igr_nlo, xpt, ratio);
          igr_nlo++;
        }
      }

      gr_werner->Set(igr_nlo);
      style(gr_werner, imu<2?2-imu:imu+1, imu<2?2-imu:imu+1);
      gr_werner->Draw("C");
      leg1->AddEntry(gr_werner, Form("#scale[0.8]{#splitline{NLO pQCD}{%s}}",leg_name[imu]), "L");
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

      style(gr_pythia, imu+2, imu<4?imu+1:imu+3);
      gr_pythia->Draw("LE");
      leg0->AddEntry(gr_pythia, Form("#scale[0.8]{#splitline{PYTHIA8}{%s}}",leg_name[imu]), "L");
    } // powheg ratio
  } // imu

  leg0->Draw();
  leg1->Draw();

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
