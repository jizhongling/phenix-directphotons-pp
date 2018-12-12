#include "QueryTree.h"
#include "DivideFunctions.h"

void draw_MissingRatio()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};
  const char *pname[3] = {"PbSc west", "PbSc east", "PbGl"};

  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root", "RECREATE");
  QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root", "RECREATE");
  QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root", "RECREATE");

  TGraphErrors *gr_miss[3];
  TGraphErrors *gr_merge[6];

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");
  THnSparse *hn_missing = (THnSparse*)f->Get("hn_missing");

  mc(0);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  mc(1, 2,1);
  legi(1, 0.2,0.7,0.4,0.9);

  for(int part=0; part<3; part++)
  {
    mcd(0);

    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    hn_missing->GetAxis(4)->SetRange(3,3);
    TH1 *h_2photon = hn_missing->Projection(1);
    hn_missing->GetAxis(3)->SetRange(2,2);
    hn_missing->GetAxis(4)->SetRange(2,2);
    TH1 *h_1photon = hn_missing->Projection(1);

    qt_miss->Fill(h_1photon, h_2photon, part);
    gr_miss[part] = DivideHisto(h_1photon, h_2photon);
    gr_miss[part]->SetNameTitle(Form("gr_%d",part), "Missing Ratio");
    aset(gr_miss[part], "p_{T}^{1#gamma} [GeV]","R", 5.,30., 0.,1.5);
    style(gr_miss[part], 20+part, 1+part);
    if(part==0)
      gr_miss[part]->Draw("AP");
    else
      gr_miss[part]->Draw("P");
    leg0->AddEntry(gr_miss[part], pname[part], "P");
    leg0->Draw();

    delete h_2photon;
    delete h_1photon;

    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    hn_missing->GetAxis(4)->SetRange(3,3);
    TH1 *h_separated = hn_missing->Projection(0);
    hn_missing->GetAxis(4)->SetRange(2,2);
    TH1 *h_merged = hn_missing->Projection(0);

    mcd(1, 1);
    gPad->SetLogy();

    qt_merge2->Fill(h_merged, h_separated, part, 0.5, 1.);
    gr_merge[part] = DivideHisto(h_merged, h_separated, 0.5, 1.);
    gr_merge[part]->SetNameTitle(Form("gr_%d",part), "Merging Ratio");
    aset(gr_merge[part], "p_{T}^{2#gamma} [GeV]","Merging ratio", 5.,30., 1e-4,10.);
    style(gr_merge[part], 20+part, 1+part);
    if(part==0)
      gr_merge[part]->Draw("AP");
    else
      gr_merge[part]->Draw("P");
    leg1->AddEntry(gr_merge[part], pname[part], "P");
    leg1->Draw();

    delete h_separated;
    delete h_merged;

    mcd(1, 2);
    gPad->SetLogy();

    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    hn_missing->GetAxis(4)->SetRange(3,3);
    h_separated = hn_missing->Projection(1);
    hn_missing->GetAxis(4)->SetRange(2,2);
    h_merged = hn_missing->Projection(1);

    qt_merge1->Fill(h_merged, h_separated, part, 0.5, 1.);
    gr_merge[part+3] = DivideHisto(h_merged, h_separated);
    gr_merge[part+3]->SetNameTitle(Form("gr_%d",part+3), "Converted Merging Ratio");
    aset(gr_merge[part+3], "p_{T}^{1#gamma} [GeV]","Converted merging ratio", 5.,30., 1e-4,10.);
    style(gr_merge[part+3], 20+part, 1+part);
    if(part==0)
      gr_merge[part+3]->Draw("AP");
    else
      gr_merge[part+3]->Draw("P");

    delete h_separated;
    delete h_merged;
  }

  c0->Print("plots/MissingRatio.pdf");
  c1->Print("plots/Merge-photon.pdf");
  qt_miss->Write();
  qt_merge1->Write();
  qt_merge2->Write();
}
