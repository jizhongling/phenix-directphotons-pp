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


void draw_ToF_Calib_ByRun() {

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
    sprintf(buf, "/phenix/plhf/zji/taxi/Run13pp510ERT/10853/data/DirectPhotonPP-%d.root", runnumber[ir]);
    TFile* f = new TFile(buf);
    //TH3* h3_tof = (TH3*)f->Get("h3_tof_raw");
    TH3* h3_tof = (TH3*)f->Get("h3_tof");

    for(Int_t is=0; is<sectors; is++) {
      TH1* h_tof = (TH1*)h3_tof->ProjectionZ("h_tof",is+1,is+1);
      h_tof->GetXaxis()->SetRangeUser(-20.,20.);
      Double_t max = h_tof->GetMaximum();
      Int_t bin1 = h_tof->FindFirstBinAbove(max/2);
      Int_t bin2 = h_tof->FindLastBinAbove(max/2);
      Double_t tof_peak = h_tof->GetBinCenter((bin1+bin2)/2);
      Double_t fwhm = h_tof->GetBinCenter(bin2) - h_tof->GetBinCenter(bin1);
      //cout << "Max=" << max << " Left bin=" << bin1 << " Right bin=" << bin2 << " Tof peak=" << tof_peak << " FWHM=" << fwhm << endl;
      g_tof_runnumber[is]->SetPoint(ir, (Double_t)runnumber[ir], tof_peak);
      g_tof_runnumber[is]->SetPointError(ir, 0., fwhm/2);
      delete h_tof;
    }

    delete h3_tof;
    delete f;

  }

  mc(0, 2,4);

  for(Int_t is=0; is<sectors; is++) {
    mcd(0, is+1);
    char buf[100];
    sprintf(buf, "Sector %d", is);
    aset(g_tof_runnumber[is]);
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

  //c0->Print("plots/ToF_NoCalib_ByRun.pdf");
  c0->Print("plots/ToF_Calib_ByRun.pdf");

}
