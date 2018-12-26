#include "GlobalVars.h"
#include "QueryTree.h"

void draw_BBCEff_Pion_Pileup()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  QueryTree *qt_pile = new QueryTree("data/BBCEff-pion-pileup.root", "RECREATE");

  QueryTree *qt_bbc = new QueryTree("data/BBCEff-pion.root");

  TF1 *fn_pol1 = new TF1("fn_pol1", "pol1");

  mc(0, 6,5);
  mc(1, 6,5);
  mc(2, 2,1);

  for(int part=0; part<2; part++)
    for(int ipt=0; ipt<npT; ipt++)
    {
      int ig = ipt*2+part;
      mcd(part, ipt+1);
      TGraphAsymmErrors *gr = qt_bbc->GraphAsymm(ig);
      mg[ig]->Draw("AP");
      mg[ig]->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      mg[ig]->Fit(fn_pol1, "Q");

      double scale = sqrt( fn_pol1->GetChisquare() / fn_pol1->GetNDF() );
      double xx = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double yy = fn_pol1->GetParameter(0);
      double eyy = fn_pol1->GetParError(0) * scale;
      if( TMath::Finite(yy+eyy) )
        qt_pile->Fill(ipt, part, xpt, yy, eyy);
    }

  for(int part=0; part<2; part++)
  {
    mcd(2, part+1);
    TGraphErrors *gr = qt_pile->Graph(part);
    gr->Set(igp[part]);
    aset(gr, "pT [GeV]","Eff", 3.,20., 0.,1.);
    style(gr, part+20, part+1);
    gr->Draw("AP");
    gr->Fit("pol0", "Q","", 3.,20.);

    gPad->Update();
    TPaveStats *st = (TPaveStats*)gr->FindObject("stats");
    st->SetY1NDC(0.6);
    st->SetY2NDC(0.8);
  }

  qt_pile->Write();
  mcw(0, "PbSc");
  mcw(1, "PbGl");
  qt_pile->Close();

  c2->Print("plots/BBCEff-pion-pileup.pdf");
}
