#include "GlobalVars.h"
#include "QueryTree.h"

void draw_Acceptance_Photon()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  QueryTree *qt_acc = new QueryTree("data/Acceptance-photon.root", "RECREATE");

  TGraphAsymmErrors *gr_acc[3];
  for(int part=0; part<3; part++)
  {
    gr_acc[part] = new TGraphAsymmErrors(npT);
    gr_acc[part]->SetName(Form("gr_%d",part));
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");
  TH1 *h_total = (TH1*)f->Get("h_photon");
  THnSparse *hn_photon = (THnSparse*)f->Get("hn_photon");

  for(int part=0; part<3; part++)
  {
    hn_photon->GetAxis(2)->SetRange(secl[part],sech[part]);
    TH1 *h_passed = hn_photon->Projection(1);
    qt_acc->Fill(h_passed, h_total, part);
    gr_acc[part]->Divide(h_passed, h_total, "n");
    delete h_passed;
  }

  mc();
  mcd();
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(int part=0; part<3; part++)
  {
    gr_acc[part]->SetTitle("Photon acceptance");
    aset(gr_acc[part], "p_{T} [GeV]","Acceptance", 4.,30., 0.,0.12);
    style(gr_acc[part], part+20, part+1);
    if(part==0)
      gr_acc[part]->Draw("AP");
    else
      gr_acc[part]->Draw("P");
    leg0->AddEntry(gr_acc[part], pname[part], "P");
  }
  leg0->Draw();

  c0->Print("plots/Acceptance-photon.pdf");
  qt_acc->Save();
}
