#include "GlobalVars.h"
#include "FitMinv.h"

void draw_ProbEff()
{
  TFile *f = new TFile("ProbEff.root");

  TH1 *h_pass[2][25];  // h_pass[part][ipt]
  TH1 *h_total[2][25];  // h_total[part][ipt]
  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<25; ipt++)
    {
      char hname[100];
      sprintf(hname,"h_pass_part%d_pt%d",part,ipt);
      h_pass[part][ipt] = (TH1*)f->Get(hname);
      sprintf(hname,"h_total_part%d_pt%d",part,ipt);
      h_total[part][ipt] = (TH1*)f->Get(hname);
    }

  TGraphErrors *gr_prob[2];
  Int_t igp[2] = {};
  for(Int_t ig=0; ig<2; ig++)
    gr_prob[ig] = new TGraphErrors(npT);

  for(Int_t part=0; part<2; part++)
    for(Int_t ic=0; ic<2; ic++)
      mc(part*2+ic, 5,5);

  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<25; ipt++)
    {
      mcd(part*2, ipt+1);
      Double_t np, enp;
      h_pass[part][ipt]->Rebin(10);
      h_pass[part][ipt]->Scale(0.5);
      h_pass[part][ipt]->SetTitle( Form("p_{T}: %4.2f-%4.2f", pTbin[ipt], pTbin[ipt+1]) );
      FitMinv(h_pass[part][ipt], np, enp);

      mcd(part*2+1, ipt+1);
      Double_t nt, ent;
      h_total[part][ipt]->Rebin(10);
      h_total[part][ipt]->Scale(0.5);
      h_total[part][ipt]->SetTitle( Form("p_{T}: %4.2f-%4.2f", pTbin[ipt], pTbin[ipt+1]) );
      FitMinv(h_total[part][ipt], nt, ent);

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
  for(Int_t part=0; part<2; part++)
    for(Int_t ic=0; ic<2; ic++)
      gROOT->ProcessLine( Form("c%d->Print(\"ProbEff-part%d-%s.pdf\")", part*2+ic, part, cond[ic]) );
  c5->Print("ProbEff.pdf");
}
