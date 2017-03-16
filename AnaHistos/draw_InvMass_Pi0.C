#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <cstdio>

using namespace std;

void draw_InvMass_Pi0()
{
  TFile *f = new TFile("/phenix/plhf/zji/taxi/Run13pp510ERT/8511/data/total.root");
  TH2F *h2_inv_mass_2photon_pi0calib_sector_raw = (TH2F*)f->Get("inv_mass_2photon_pi0calib_sector_raw");
  TH2F *h2_inv_mass_2photon_pi0calib_sector = (TH2F*)f->Get("inv_mass_2photon_pi0calib_sector");

  TCanvas *c = new TCanvas("c", "Canvas", 1200, 600);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c->Divide(2,1);

  c->cd(1);
  h2_inv_mass_2photon_pi0calib_sector_raw->SetTitle("m_{inv} before Calibration");
  h2_inv_mass_2photon_pi0calib_sector_raw->Draw("colz");

  c->cd(2);
  h2_inv_mass_2photon_pi0calib_sector->SetTitle("m_{inv} after Calibration");
  h2_inv_mass_2photon_pi0calib_sector->Draw("colz");

  c->Print("InvMass_Pi0-1.pdf");
  c->Clear("D");

  c->cd(1);
  TH1D *hp_inv_mass_2photon_pi0calib_sector_raw = (TH1D*)h2_inv_mass_2photon_pi0calib_sector_raw->ProjectionX("_px",1,1);
  hp_inv_mass_2photon_pi0calib_sector_raw->SetTitle("m_{inv} mass before Calibration for Sector 0");
  hp_inv_mass_2photon_pi0calib_sector_raw->Fit("gaus", "", "", 0.12, 0.15);
  
  c->cd(2);
  TH1D *hp_inv_mass_2photon_pi0calib_sector = (TH1D*)h2_inv_mass_2photon_pi0calib_sector->ProjectionX("_px",1,1);
  hp_inv_mass_2photon_pi0calib_sector->SetTitle("m_{inv} mass after Calibration for Sector 0");
  hp_inv_mass_2photon_pi0calib_sector->Fit("gaus", "", "", 0.12, 0.15);
  
  c->Print("InvMass_Pi0-2.pdf");
  c->Clear("D");

  TH2F *h2_tof_sector_raw = (TH2F*)f->Get("tof_sector_raw");
  TH2F *h2_tof_sector = (TH2F*)f->Get("tof_sector");

  c->cd(1);
  h2_tof_sector_raw->SetTitle("ToF before Calibration");
  h2_tof_sector_raw->Draw("colz");

  c->cd(2);
  h2_tof_sector->SetTitle("ToF after Calibration");
  h2_tof_sector->Draw("colz");

  c->Print("ToF_Calib_Sector-1.pdf");
  c->Clear("D");

  c->cd(1);
  TH1D *hp_tof_sector_raw = (TH1D*)h2_tof_sector_raw->ProjectionX("_px",1,1);
  hp_tof_sector_raw->SetTitle("ToF before Calibration for Sector 0");
  hp_tof_sector_raw->Fit("gaus", "", "", -30., 10.);

  c->cd(2);
  TH1D *hp_tof_sector = (TH1D*)h2_tof_sector->ProjectionX("_px",1,1);
  hp_tof_sector->SetTitle("ToF after Calibration for Sector 0");
  hp_tof_sector->Fit("gaus", "", "", -20., 20.);

  c->Print("ToF_Calib_Sector-2.pdf");
  c->Clear("D");
}
