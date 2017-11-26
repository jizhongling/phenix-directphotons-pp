#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TH1D.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TCanvas.h>
#include <TStyle.h>

using namespace std;


void draw_ToF_Calib_Tower() {

  ofstream fout("ToF_Calib_Tower.txt");
  TFile* f = new TFile("/phenix/plhf/zji/taxi/Run13pp510MinBias/8328/data/total.root");
  TH2I* h_hitmap_tof_energy_PbSc = (TH2I*)f->Get("hitmap_tof_energy_PbSc");
  TH2I* h_hitmap_tof_energy_PbGl = (TH2I*)f->Get("hitmap_tof_energy_PbGl");
  TH2I* h_hitmap_tof_tower = (TH2I*)f->Get("hitmap_tof_tower");

  TCanvas* c = new TCanvas("c", "Canvas", 1800, 600);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  c->Divide(3,1);

  c->cd(1);
  h_hitmap_tof_energy_PbSc->Draw("colz");

  c->cd(2);
  h_hitmap_tof_energy_PbGl->Draw("colz");

  c->cd(3);
  fout << "towerid    ToF peak" << endl; 
  for(Int_t id=0; id<25000; id++) {
    TH1D* hp_hitmap_tof_tower = (TH1D*)h_hitmap_tof_tower->ProjectionY("_py",id+1,id+1);
    TSpectrum* peak = new TSpectrum();
    Int_t nfound = peak->Search(hp_hitmap_tof_tower,2.,"nodraw");
    if(nfound) {
      Float_t peakX = *peak->GetPositionX();
      hp_hitmap_tof_tower->Fit("gaus","Q0","",peakX-3.,peakX+3.);
      TF1* f_hitmap_tof_tower = hp_hitmap_tof_tower->GetFunction("gaus");
      Double_t mean_hitmap_tof_tower = f_hitmap_tof_tower->GetParameter(1);
      fout << id << "    " << mean_hitmap_tof_tower << endl;
      f_hitmap_tof_tower->Delete();
    }
    else {
      fout << id << "    " << "-100." << endl;
    }
    delete peak;
    hp_hitmap_tof_tower->Delete();
  }
  fout.close();

  TH1D* hp_hitmap_tof_tower = (TH1D*)h_hitmap_tof_tower->ProjectionY("_py",100,100);
  TSpectrum* peak = new TSpectrum();
  peak->Search(hp_hitmap_tof_tower,2.,"nodraw");
  Float_t peakX = *peak->GetPositionX();
  hp_hitmap_tof_tower->Fit("gaus","","",peakX-3.,peakX+3.);
  TF1* f_hitmap_tof_tower = hp_hitmap_tof_tower->GetFunction("gaus");
  Double_t mean_hitmap_tof_tower = f_hitmap_tof_tower->GetParameter(1);
  cout << "peakX=" << peakX << endl;
  cout << "ToF mean=" << mean_hitmap_tof_tower << endl;

  c->Print("plots/ToF_Calib_Tower.pdf");

}
