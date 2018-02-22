#include "GlobalVars.h"

void draw_BBCEff_Photon_Pileup()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const Int_t secl[2] = {1, 7};
  const Int_t sech[2] = {6, 8};

  TMultiGraph *mg[npT*2];
  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      Int_t ig = ipt*2+part;
      mg[ig] = new TMultiGraph();
    }

  for(Int_t i=0; i<1; i++)
  {
    TFile *f = new TFile(Form("pileup/BBCEff-photon-%d.root",i));
    if( f->IsZombie() ) continue;

    for(Int_t part=0; part<2; part++)
      for(Int_t ipt=0; ipt<npT; ipt++)
      {
        Int_t ig = ipt*2+part;
        TGraphAsymmErrors *gr = (TGraphAsymmErrors*)f->Get(Form("gr_%d",ig));
        mg[ig]->Add(gr);
      }
  }

  TGraphErrors *gr_trig[2];
  Int_t igp[2] = {};
  for(Int_t part=0; part<2; part++)
    gr_trig[part] = new TGraphErrors(npT);

  TF1 *fn_pol1 = new TF1("fn_pol1", "pol1");

  mc(0, 6,5);
  mc(1, 6,5);
  mc(2, 2,1);

  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      Int_t ig = ipt*2+part;
      mcd(part, ipt+1);
      mg[ig]->Draw("AP");
      mg[ig]->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      mg[ig]->Fit(fn_pol1, "Q");

      Double_t scale = sqrt( fn_pol1->GetChisquare() / fn_pol1->GetNDF() );
      Double_t xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      Double_t yy = fn_pol1->GetParameter(0);
      Double_t eyy = fn_pol1->GetParError(0) * scale;
      if( yy > 0. && eyy > 0. && eyy < TMath::Infinity() )
      {
        gr_trig[part]->SetPoint(igp[part], xx, yy);
        gr_trig[part]->SetPointError(igp[part], 0., eyy);
        igp[part]++;
      }
    }

  for(Int_t part=0; part<2; part++)
  {
    mcd(2, part+1);
    gr_trig[part]->Set(igp[part]);
    aset(gr_trig[part], "pT [GeV]","Eff", 3.,20., 0.,1.);
    style(gr_trig[part], part+20, part+1);
    gr_trig[part]->Draw("AP");
    gr_trig[part]->Fit("pol0", "Q","", 3.,20.);

    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr_trig[part]->FindObject("stats");
    st->SetY1NDC(0.6);
    st->SetY2NDC(0.8);
  }

  TFile *f_out = new TFile("data/BBCEff-photon-pileup.root", "RECREATE");
  mcw(0, "PbSc");
  mcw(1, "PbGl");
  f_out->Close();

  c2->Print("plots/BBCEff-photon-pileup.pdf");
}
