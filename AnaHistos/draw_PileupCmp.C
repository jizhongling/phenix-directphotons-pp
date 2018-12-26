#include "Pileup.h"

void draw_PileupCmp()
{
  int gn_mine[2] = {}, gn_sasha[2] = {};
  double runno_mine[2][1000] = {}, runno_sasha[2][1000] = {};
  double npi0_mine[2][1000] = {}, npi0_sasha[2][1000] = {};

  for(int i=0; i<44; i++)
  {
    TFile *f_mine = new TFile(Form("histos/Pileup-%d.root",i));
    if( f_mine->IsZombie() ) continue;

    for(int igr=0; igr<2; igr++)
    {
      TGraphErrors *gr = (TGraphErrors*)f_mine->Get(Form("gr_run_%d",2*8*npT+4+igr));
      for(int ip=0; ip<gr->GetN(); ip++)
      {
        gr->GetPoint(ip, runno_mine[igr][gn_mine[igr]], npi0_mine[igr][gn_mine[igr]]);
        gn_mine[igr]++;
      }
    }

    delete f_mine;
  }

  TFile *f_sasha = new TFile("data/Pileup-Sasha-CVS.root");
  TFile *f_sasha0 = new TFile("data/Pileup-Sasha-TAXI.root");

  for(int igr=0; igr<2; igr++)
  {
    TGraphErrors *gr = (TGraphErrors*)f_sasha->Get(Form("gr_run_%d",2+igr));
    TGraphErrors *gr0 = (TGraphErrors*)f_sasha0->Get(Form("gr_run_%d",2+igr));
    for(int ip=0; ip<gr->GetN(); ip++)
    {
      gr->GetPoint(ip, runno_sasha[igr][gn_sasha[igr]], npi0_sasha[igr][gn_sasha[igr]]);
      gn_sasha[igr]++;
      //gr0->GetPoint(ip, runno_mine[igr][gn_mine[igr]], npi0_mine[igr][gn_mine[igr]]);
      //gn_mine[igr]++;
    }
  }

  delete f_sasha;
  delete f_sasha0;

  for(int igr=0; igr<2; igr++)
  {
    for(int j=0; j<gn_sasha[igr]; j++)
    {
      bool matched = false;
      for(int i=0; i<gn_mine[igr]; i++)
        if( TMath::Abs(runno_mine[igr][i]-runno_sasha[igr][j]) < 0.1 )
          matched = true;
      if(!matched)
        cout << "Part " << igr << " Runnumber = " << runno_sasha[igr][j] << " not matched!!!" << endl;
    }
  }

  mc();
  mcd();

  TGraph *gr_ratio[2];
  for(int igr=0; igr<2; igr++)
  {
    gr_ratio[igr] = new TGraph(1000);
    int igp1 = 0;
    for(int i=0; i<gn_mine[igr]; i++)
      for(int j=0; j<gn_sasha[igr]; j++)
        if( TMath::Abs(runno_mine[igr][i]-runno_sasha[igr][j]) < 0.1 )
        {
          double xx = runno_mine[igr][i];
          double yy = npi0_mine[igr][i] / npi0_sasha[igr][j];
          gr_ratio[igr]->SetPoint(igp1, xx, yy);
          igp1++;
        }
    gr_ratio[igr]->Set(igp1);
    aset(gr_ratio[igr], "Runnumber","Mine/Sasha", 386000.,400000., 0.9, 1.1);
    style(gr_ratio[igr], igr+24, igr+1);
    cout << "gn_mine[" << igr << "] = " << gn_mine[igr] << endl
      << "gn_sasha[" << igr << "] = " << gn_sasha[igr] << endl
      << "gn_ratio[" << igr << "] = " << gr_ratio[igr]->GetN() << endl;
  }
  gr_ratio[0]->Draw("AP");
  gr_ratio[1]->Draw("P");

  c0->Print("plots/PileupCmp.pdf");
}
