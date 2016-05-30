#include <iostream>
#include <fstream>
#include <stdio.h>

#include <TFile.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>

using namespace std;


void draw_ToF_Calib_Sector() {

  const Int_t sectors = 8;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");
  Int_t nrun = 0;
  Int_t runnumber[1024];
  while(fin >> runnumber[nrun]) nrun++;
  fin.close();

  TGraphErrors* g_tof_runnumber[sectors];
  for(Int_t is=0; is<sectors; is++)
    g_tof_runnumber[is] = new TGraphErrors(nrun);

  for(Int_t ir=0; ir<nrun; ir++) {

    char buf[100];
    sprintf(buf, "/phenix/plhf/zji/taxi/Run13pp510ERT/8511/data/DirectPhotonPP-%d.root", runnumber[ir]);
    TFile* f = new TFile(buf);
    //TH2I* h_hitmap_tof_sector = (TH2I*)f->Get("tof_sector_raw");
    TH2I* h_hitmap_tof_sector = (TH2I*)f->Get("tof_sector");

    for(Int_t is=0; is<sectors; is++) {
      TH1D* hp_hitmap_tof_sector = (TH1D*)h_hitmap_tof_sector->ProjectionX("_px",is+1,is+1);
      hp_hitmap_tof_sector->GetXaxis()->SetRangeUser(-20.,20.);
      Double_t max = hp_hitmap_tof_sector->GetMaximum();
      Int_t bin1 = hp_hitmap_tof_sector->FindFirstBinAbove(max/2);
      Int_t bin2 = hp_hitmap_tof_sector->FindLastBinAbove(max/2);
      Double_t tof_peak = hp_hitmap_tof_sector->GetBinCenter((bin1+bin2)/2);
      Double_t fwhm = hp_hitmap_tof_sector->GetBinCenter(bin2) - hp_hitmap_tof_sector->GetBinCenter(bin1);
      //cout << "Max=" << max << " Left bin=" << bin1 << " Right bin=" << bin2 << " Tof peak=" << tof_peak << " FWHM=" << fwhm << endl;
      g_tof_runnumber[is]->SetPoint(ir, (Double_t)runnumber[ir], tof_peak);
      g_tof_runnumber[is]->SetPointError(ir, 0., fwhm/2);
      hp_hitmap_tof_sector->Delete();
    }

    h_hitmap_tof_sector->Delete();
    delete f;

  }

  TCanvas* c = new TCanvas("c", "Canvas", 1200, 2400);
  gStyle->SetOptStat(1);
  c->Divide(2,4);

  for(Int_t is=0; is<sectors; is++) {
    c->cd(is+1);
    char buf[100];
    sprintf(buf, "Sector %d", is);
    g_tof_runnumber[is]->SetTitle(buf);
    g_tof_runnumber[is]->GetXaxis()->SetTitle("Runnumber");
    g_tof_runnumber[is]->GetXaxis()->SetLimits(387000.,399000.);
    g_tof_runnumber[is]->GetYaxis()->SetTitle("ToF peak [ns]");
    g_tof_runnumber[is]->GetYaxis()->SetRangeUser(-50.,50.);
    g_tof_runnumber[is]->SetMarkerColor(1);
    g_tof_runnumber[is]->SetMarkerSize(1);
    g_tof_runnumber[is]->SetMarkerStyle(21);
    g_tof_runnumber[is]->Draw("AP");
  }

  //c->Print("ToF_NoCalib_Sector.pdf");
  c->Print("ToF_Calib_Sector.pdf");

}
