#include "GlobalVars.h"
#include "QueryTree.h"
#include "DivideFunctions.h"

void draw_Iso2Inc(const int pwhg = 0, const int ipwhg = 0)
{
  const int sector = 3;  // PbSc west: 0; PbSc east: 1; PbGl: 2; Combined: 3
  const char *pwhg_type[4] = {"with MPI", "QED-QCD veto", "without MPI", "pure hard"};
  const char *suffix[4] = {"-pwhg", "-qedqcd", "-nompi", "-purehard"};
  const char *jetphox_fname[3] = {"halfpt", "onept", "twopt"};
  const char *jetphox_scale[3] = {"p_{T}/2", "p_{T}", "2p_{T}"};

  QueryTree *qt_pion = new QueryTree("data/CrossSection-pion.root");
  QueryTree *qt_photon = new QueryTree("data/CrossSection-photon.root");
  QueryTree *qt_isophoton = new QueryTree("data/CrossSection-isophoton.root");
  QueryTree *qt_jetphox = new QueryTree("data/JetphoxRatio.root");
  QueryTree *qt_sys = new QueryTree("data/CrossSection-syserr.root");

  if(pwhg == 0)
  {
    TFile *f_pythia = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-photon.root");
    THnSparse *hn_hadron = (THnSparse*)f_pythia->Get("hn_hadron");
    hn_hadron->GetAxis(3)->SetRange(1,1);  // prompt photons
    hn_hadron->GetAxis(2)->SetRange(1,1);  // |eta| < 0.25
    hn_hadron->GetAxis(5)->SetRange(1,2);  // inclusive
    TH1 *h_photon = hn_hadron->Projection(0);
    h_photon->SetName("h_photon_eta025");
    hn_hadron->GetAxis(5)->SetRange(2,2);  // isolated
    TH1 *h_isophoton = hn_hadron->Projection(0);
    h_isophoton->SetName("h_isophoton_eta025");
  }
  else if(pwhg == 1)
  {
    TFile *f_pythia = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaPowheg-histo%s.root",ipwhg?suffix[ipwhg]:""));
    TH1 *h_photon = (TH1*)f_pythia->Get("hard0_iso0_rap0_id0");
    TH1 *h_isophoton = (TH1*)f_pythia->Get("hard0_iso1_rap0_id0");
  }
  else
  {
    cout << "Wrong input" << endl;
    return;
  }

  TGraphErrors *gr[2];
  int igp[2] = {};
  for(int iph=0; iph<2; iph++)
    gr[iph] = new TGraphErrors(npT);
  TGraphErrors *gr_sys = new TGraphErrors(npT);

  for(int ipt=12; ipt<npT; ipt++)
  {
    double xpt, incl, eincl, iso, eiso;

    if( qt_pion->Query(ipt, sector, xpt, incl, eincl) &&
        qt_pion->Query(ipt, sector+4, xpt, iso, eiso) )
    {
      double yy = iso / incl;
      double eyy = yy * sqrt( pow(eincl/incl,2) + pow(eiso/iso,2) );
      if( TMath::Finite(yy+eyy) )
      {
        gr[0]->SetPoint(igp[0], xpt, yy);
        gr[0]->SetPointError(igp[0], 0., eyy);
        igp[0]++;
      }
    }

    if( qt_photon->Query(ipt, sector, xpt, incl, eincl) &&
        qt_isophoton->Query(ipt, sector, xpt, iso, eiso) )
    {
      double yy = iso / incl;
      double eyy = yy * sqrt( pow(eincl/incl,2) + pow(eiso/iso,2) );
      if( TMath::Finite(yy+eyy) )
      {
        gr[1]->SetPoint(igp[1], xpt, yy);
        gr[1]->SetPointError(igp[1], 0., eyy);
        double xsec[2], sys[2];
        for(int i=0; i<2; i++)
          qt_sys->Query(ipt, i, xpt, xsec[i], sys[i]);
        double yy = xsec[1] / xsec[0];
        double eyy = yy * sqrt( pow(sys[0]/xsec[0],2) + pow(sys[1]/xsec[1],2) );
        gr_sys->SetPoint(igp[1], xpt, yy);
        gr_sys->SetPointError(igp[1], 0., eyy);
        igp[1]++;
      }
    }
  }

  mc();
  mcd();
  legi(0, 0.22,0.8,0.85,0.9);
  leg0->SetNColumns(3);
  legi(1, 0.22,0.65,0.42,0.82);
  for(int iph=0; iph<2; iph++)
  {
    gr[iph]->Set(igp[iph]);
    aset(gr[iph], "p_{T} [GeV/c]", "#frac{Isolated}{Inclusive}", 4.9,30.1, 0.,1.6);
    style(gr[iph], iph+20, iph+1);
    if(iph==0)
    {
      gr[iph]->SetTitle("");
      gr[iph]->Draw("AP");
    }
    else
    {
      gr[iph]->Draw("P");
      style(gr_sys, 1, iph+1);
      gr_sys->SetLineWidth(2);
      gr_sys->Draw("[]");
    }
  }
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
  style(gr_pythia, 1, 4);
  gr_pythia->Draw("LE");
  leg0->AddEntry(gr[0], "#scale[0.9]{Data #pi^{0}}", "P");
  leg0->AddEntry(gr[1], "#scale[0.9]{Data #gamma_{dir}}", "P");
  leg0->AddEntry(gr_pythia, Form("#scale[0.7]{#splitline{Pythia #gamma_{dir}}{%s}}",pwhg_type[ipwhg]), "L");
  for(int imu=0; imu<3; imu++)
  {
    TFile *f_iso = new TFile( Form("data/isoprompt-x2000-ct14-%s.root",jetphox_fname[imu]) );
    TFile *f_inc = new TFile( Form("data/incprompt-x2000-ct14-%s.root",jetphox_fname[imu]) );
    TH1 *h_iso = (TH1*)f_iso->Get("hp41");
    TH1 *h_inc = (TH1*)f_inc->Get("hp41");
    TGraphErrors *gr_jetphox = new TGraphErrors(npT);
    int igr_jetphox = 0;
    for(int ipt=12; ipt<npT; ipt++)
    {
      double xpt = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      int bin_th = h_iso->GetXaxis()->FindBin(xpt);
      double ratio = h_iso->GetBinContent(bin_th) / h_inc->GetBinContent(bin_th);
      double eratio = ratio * sqrt( pow(h_iso->GetBinError(bin_th)/h_iso->GetBinContent(bin_th),2) + pow(h_inc->GetBinError(bin_th)/h_inc->GetBinContent(bin_th),2) );
      if( TMath::Finite(ratio+eratio) )
      {
        gr_jetphox->SetPoint(igr_jetphox, xpt, ratio);
        gr_jetphox->SetPointError(igr_jetphox, 0., eratio);
        igr_jetphox++;
      }
    }
    gr_jetphox->Set(igr_jetphox);
    delete f_iso;
    delete f_inc;
    style(gr_jetphox, 2, imu+1);
    gr_jetphox->Draw("LE");
    leg1->AddEntry(gr_jetphox, Form("#scale[0.7]{JETPHOX #gamma_{dir} %s}",jetphox_scale[imu]), "L");
  }
  leg0->Draw();
  leg1->Draw();

  const char *outfile = Form("plots/Iso2Inc%s", pwhg?suffix[ipwhg]:"-pythia6");
  c0->Print(Form("%s.pdf", outfile));
  if(pwhg == 1 && ipwhg != 1)
  {
    char *cmd = Form("preliminary.pl --input=%s.pdf --output=%s-prelim.pdf --x=320 --y=200 --scale=0.8", outfile,outfile);
    system(cmd);
    if(ipwhg != 0)
    {
      cmd = Form("rm %s.pdf", outfile);
      system(cmd);
    }
  }
}
