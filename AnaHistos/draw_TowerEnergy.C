#include "GlobalVars.h"

void draw_TowerEnergy()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TFile *f = new TFile("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/13173/data/PhotonHistos-histo.root");
  THnSparse *hn_etwr = (THnSparse*)f->Get("hn_etwr");

  hn_etwr->GetAxis(7)->SetRange(1,1);
  hn_etwr->GetAxis(1)->SetRange(26,-1);
  hn_etwr->GetAxis(2)->SetRange(4,4);
  hn_etwr->GetAxis(3)->SetRange(4,4);
  for(int sec=0; sec<8; sec++)
  {
    hn_etwr->GetAxis(0)->SetRange(sec+1,sec+1);
    TH1 *h_cut = hn_etwr->Projection(6);
    double ncount = 0.;
    for(int ic=2; ic<=8; ic+=2)
      ncount += h_cut->GetBinContent(ic);
    cout << "Sector " << sec << ": " << ncount << endl;
    delete h_cut;
  }
  return;

  for(int ic=0; ic<2; ic++)
    for(int part=0; part<3; part++)
    {
      mc(ic*3+part, 6,5);
      mc(ic*3+part+6, 6,5);
      for(int ipt=0; ipt<npT; ipt++)
      {
        hn_etwr->GetAxis(7)->SetRange(ic+1,ic+1);
        hn_etwr->GetAxis(0)->SetRange(secl[part],sech[part]);
        hn_etwr->GetAxis(1)->SetRange(ipt+1,ipt+1);

        mcd(ic*3+part, ipt+1);
        gPad->SetLogz();
        TH2 *h2_yz = hn_etwr->Projection(3,2);
        h2_yz->SetTitle( Form("p_{T}: %.1f-%.1f",pTbin[ipt],pTbin[ipt+1]) );
        h2_yz->DrawCopy("COLZ");
        delete h2_yz;

        mcd(ic*3+part+6, ipt+1);
        TH1 *h_tof = hn_etwr->Projection(4);
        h_tof->SetTitle( Form("p_{T}: %.1f-%.1f",pTbin[ipt],pTbin[ipt+1]) );
        aset(h_tof);
        style(h_tof, 20, 1);
        h_tof->DrawCopy();
        delete h_tof;
      }
    }

  TFile *f_out = new TFile("data/ETower.root", "RECREATE");
  for(int ic=0; ic<2; ic++)
    for(int part=0; part<3; part++)
    {
      mcw( ic*3+part, Form("part%d-bbc%d-yz",part,ic) );
      mcw( ic*3+part+6, Form("part%d-bbc%d-tof",part,ic) );
    }
  f_out->Close();
}
