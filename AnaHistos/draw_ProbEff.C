#include "GlobalVars.h"
#include "FitMinv.h"

void draw_ProbEff()
{
  TFile *f = new TFile("Pi0PP-histo.root");

  TH1 *h_pass[3][25];  // h_pass[is][ipt]
  TH1 *h_total[3][25];  // h_total[is][ipt]
  for(Int_t is=0; is<3; is++)
    for(Int_t ipt=0; ipt<25; ipt++)
    {
      char hname[100];
      sprintf(hname, "mchist_s%d_pt%02d_tp", is, ipt);
      h_pass[is][ipt] = (TH1*)f->Get(hname);
      sprintf(hname, "mchist_s%d_pt%02d_t", is, ipt);
      h_total[is][ipt] = (TH1*)f->Get(hname);
    }

  for(Int_t ipt=0; ipt<25; ipt++)
  {
    h_pass[0][ipt]->Add(h_pass[1][ipt]);
    h_total[0][ipt]->Add(h_total[1][ipt]);
  }

  TGraphErrors *gr_prob[2];
  Int_t igp[2] = {};
  for(Int_t part=0; part<2; part++)
  {
    gr_prob[part] = new TGraphErrors(npT);
    gr_prob[part]->SetName(Form("gr_%d",part));
  }

  for(Int_t part=0; part<2; part++)
    for(Int_t ic=0; ic<2; ic++)
      mc(part*2+ic, 5,5);

  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<25; ipt++)
    {
      mcd(part*2, ipt+1);
      Double_t np, enp;
      h_pass[part*2][ipt]->Rebin(10);
      h_pass[part*2][ipt]->Scale(0.5);
      h_pass[part*2][ipt]->SetTitle( Form("p_{T}: %4.2f-%4.2f", pTbin[ipt], pTbin[ipt+1]) );
      FitMinv(h_pass[part*2][ipt], np, enp);

      mcd(part*2+1, ipt+1);
      Double_t nt, ent;
      h_total[part*2][ipt]->Rebin(10);
      h_total[part*2][ipt]->Scale(0.5);
      h_total[part*2][ipt]->SetTitle( Form("p_{T}: %4.2f-%4.2f", pTbin[ipt], pTbin[ipt+1]) );
      FitMinv(h_total[part*2][ipt], nt, ent);

      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t yy = np / nt;
      Double_t eyy = yy * sqrt( pow(ent/nt,2.) + pow(enp/np,2.) );
      if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_prob[part]->SetPoint(igp[part], xx, yy);
        gr_prob[part]->SetPointError(igp[part], 0., eyy);
        igp[part]++;
      }
    }

  mc(5, 2,1);

  gr_prob[0]->SetTitle("Prob Eff for PbSc");
  gr_prob[1]->SetTitle("Prob Eff for PbGl");

  for(Int_t part=0; part<2; part++)
  {
    mcd(5, part+1);
    gr_prob[part]->Set(igp[part]);
    aset(gr_prob[part], "p_{T} [GeV]","Prob Eff", 0.,20., 0.8,1.1);
    style(gr_prob[part], 24, kRed);
    gr_prob[part]->Draw("AP");
  }

  const char* cond[2] = {"pass", "total"};
  TFile *f_out = new TFile("ProbEff.root", "RECREATE");
  for(Int_t part=0; part<2; part++)
  {
    gr_prob[part]->Write();
    for(Int_t ic=0; ic<2; ic++)
      mcw(part*2+ic, Form("part%d-%s",part,cond[ic]));
  }
  f_out->Close();

  c5->Print("ProbEff.pdf");
}
