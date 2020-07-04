#include "GlobalVars.h"
#include "QueryTree.h"
#include "DivideFunctions.h"

void draw_Iso2Inc(const int pwhg = 0)
{
  const int sector = 3;  // PbSc west: 0; PbSc east: 1; PbGl: 2; Combined: 3
  const char *jetphox_fname[3] = {"p_{T}/2", "p_{T}", "2p_{T}"};

  QueryTree *qt_pion = new QueryTree("data/CrossSection-pion.root");
  QueryTree *qt_photon = new QueryTree("data/CrossSection-photon.root");
  QueryTree *qt_isophoton = new QueryTree("data/CrossSection-isophoton.root");
  QueryTree *qt_jetphox = new QueryTree("data/JetphoxRatio.root");
  QueryTree *qt_sys = new QueryTree("data/CrossSection-syserr.root");

  if(pwhg == 0)
  {
    TFile *f_pythia = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-minbias.root");
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
    TFile *f_pythia = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaPowheg-histo-x703.root");
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

  for(int ipt=0; ipt<npT; ipt++)
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
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  for(int iph=0; iph<2; iph++)
  {
    gr[iph]->Set(igp[iph]);
    aset(gr[iph], "p_{T} [GeV]", "Iso/Inc", 6.,30., 0.,1.4);
    style(gr[iph], iph+20, iph+1);
    if(iph==0)
    {
      gr[iph]->SetTitle("Isolated/Inclusive ratio");
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
  style(gr_pythia, 25, 2);
  gr_pythia->Draw("P");
  leg0->AddEntry(gr[0], "Data #pi^{0}", "P");
  leg0->AddEntry(gr[1], "Data #gamma_{dir}", "P");
  leg0->AddEntry(gr_pythia, "Pythia #gamma_{dir}", "P");
  for(int imu=0; imu<3; imu++)
  {
    TGraphErrors *gr_jetphox = qt_jetphox->Graph(imu);
    style(gr_jetphox, imu+20, imu+1);
    gr_jetphox->Draw("LE");
    leg0->AddEntry(gr_jetphox, Form("NLO %s",jetphox_fname[imu]), "L");
  }
  leg0->Draw();

  c0->Print(Form("plots/Iso2Inc%s.pdf",pwhg?"-pwhg":""));
}
