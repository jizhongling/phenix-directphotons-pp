#include "QueryTree.h"

void draw_ProbEff_PISA()
{
  const char *pname[2] = {"PbSc", "PbGl"};
  const int secl[2] = {1, 7};
  const int sech[2] = {6, 8};

  QueryTree *qt_prob = new QueryTree("data/ProbEff-PISA.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/HadronResponse-histo-photon.root");
  THnSparse *hn_prob = (THnSparse*)f->Get("hn_prob_photon");
  TAxis *axis_pt = hn_prob->GetAxis(0);
  TAxis *axis_sec = hn_prob->GetAxis(1);
  TAxis *axis_photon = hn_prob->GetAxis(2);
  TAxis *axis_prob = hn_prob->GetAxis(3);

  mc(0, 2,1);

  for(int part=0; part<2; part++)
  {
    axis_photon->SetRange(2,2);
    axis_sec->SetRange(secl[part],sech[part]);
    axis_prob->SetRange(1,2);
    TH1 *h_total = hn_prob->Projection(0);
    axis_prob->SetRange(2,2);
    TH1 *h_passed = hn_prob->Projection(0);
    qt_prob->Fill(h_passed, h_total, part);

    mcd(0, part+1);
    TGraphErrors *gr = qt_prob->Graph(part);
    gr->SetTitle( Form("Prob Eff for %s",pname[part]) );
    aset(gr, "p_{T} (GeV/c)","Prob Eff", 5.1,30., 0.95,1.01);
    style(gr, part+20, part+1);
    gr->Draw("AP");

    delete h_passed;
    delete h_total;
  }

  c0->Print("plots/ProbEff-PISA.pdf");
  qt_prob->Save();
}
