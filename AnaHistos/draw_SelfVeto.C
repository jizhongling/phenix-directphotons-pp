#include "GlobalVars.h"
#include "QueryTree.h"

void draw_SelfVeto()
{
  const int secl[4] = {1, 5, 7, 1};
  const int sech[4] = {4, 6, 8, 8};
  const char *pname[4] = {"PbSc west", "PbSc east", "PbGl", "Combined"};

  QueryTree *qt_veto = new QueryTree("data/SelfVeto.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");
  TH3 *h3_isopair = (TH3*)f->Get("h3_isoeta");

  mc();
  mcd();
  legi(0, 0.5,0.7,0.8,0.9);

  for(int part=0; part<4; part++)
  {
    int npT_rebin = part<3 ? npT : npT_pol;
    double *pTrebin = part<3 ? pTbin : pTbin_pol;

    TH1 *h_total = h3_isopair->ProjectionX("h_total0", secl[part],sech[part], 1,2);
    TH1 *h_passed = h3_isopair->ProjectionX("h_passed0", secl[part],sech[part], 2,2);
    h_total = h_total->Rebin(npT_rebin, "h_total", pTrebin);
    h_passed = h_passed->Rebin(npT_rebin, "h_passed", pTrebin);

    qt_veto->Fill(h_passed, h_total, part);
    TGraphErrors *gr = qt_veto->Graph(part);
    gr->SetName(Form("gr_%d",part));

    aset(gr, "p_{T} (GeV/c)","#frac{isoboth}{isopair}", 0.,30., 0.,0.6);
    style(gr, 20+part, 1+part);
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    leg0->AddEntry(gr, pname[part], "P");

    delete h_total;
    delete h_passed;
  }
  leg0->Draw();

  c0->Print("plots/SelfVeto.pdf");
  qt_veto->Save();
}
