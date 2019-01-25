#include "QueryTree.h"

void draw_MissingRatio()
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};
  const char *pname[3] = {"PbSc west", "PbSc east", "PbGl"};

  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root", "RECREATE");
  QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root", "RECREATE");
  QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root", "RECREATE");
  QueryTree *qt_ptratio = new QueryTree("data/PtRatio-sim.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");
  THnSparse *hn_missing = (THnSparse*)f->Get("hn_missing");

  mc(0);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  mc(1, 2,1);
  legi(1, 0.2,0.7,0.4,0.9);
  mc(2);

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
    TGraphErrors *gr_miss = qt_miss->Graph(part);
    gr_miss->SetNameTitle(Form("gr_%d",part), "Missing Ratio");
    aset(gr_miss, "p_{T}^{1#gamma} [GeV]","R", 5.,30., 0.,1.5);
    style(gr_miss, 20+part, 1+part);
    if(part==0)
      gr_miss->Draw("AP");
    else
      gr_miss->Draw("P");
    leg0->AddEntry(gr_miss, pname[part], "P");
    leg0->Draw();

    delete h_2photon;
    delete h_1photon;

    mcd(1, 1);
    gPad->SetLogy();

    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    hn_missing->GetAxis(4)->SetRange(3,3);
    h_separated = hn_missing->Projection(1);
    hn_missing->GetAxis(4)->SetRange(2,2);
    h_merged = hn_missing->Projection(1);

    qt_merge1->Fill(h_merged, h_separated, part);
    TGraphErrors *gr_merge1 = qt_merge1->Graph(part);
    gr_merge1->SetNameTitle(Form("gr_%d",part+3), "Converted Merging Ratio");
    aset(gr_merge1, "p_{T}^{1#gamma} [GeV]","Converted merging ratio", 5.,30., 1e-4,10.);
    style(gr_merge1, 20+part, 1+part);
    if(part==0)
      gr_merge1->Draw("AP");
    else
      gr_merge1->Draw("P");
    leg1->AddEntry(gr_merge1, pname[part], "P");
    leg1->Draw();

    delete h_separated;
    delete h_merged;

    mcd(1, 2);
    gPad->SetLogy();

    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    hn_missing->GetAxis(4)->SetRange(3,3);
    TH1 *h_separated = hn_missing->Projection(0);
    hn_missing->GetAxis(4)->SetRange(2,2);
    TH1 *h_merged = hn_missing->Projection(0);

    qt_merge2->Fill(h_merged, h_separated, part);
    TGraphErrors *gr_merge2 = qt_merge2->Graph(part);
    gr_merge2->SetNameTitle(Form("gr_%d",part), "Merging Ratio");
    aset(gr_merge2, "p_{T}^{2#gamma} [GeV]","Merging ratio", 5.,30., 1e-4,20.);
    style(gr_merge2, 20+part, 1+part);
    if(part==0)
      gr_merge2->Draw("AP");
    else
      gr_merge2->Draw("P");

    delete h_separated;
    delete h_merged;

    mcd(2);

    hn_missing->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_missing->GetAxis(3)->SetRange(3,3);
    //hn_missing->GetAxis(4)->SetRange(3,3);
    h_pt1photon = hn_missing->Projection(1);
    h_pt2photon = hn_missing->Projection(0);

    qt_ptratio->Fill(h_pt2photon, h_pt1photon, part);
    TGraphErrors *gr_ptratio = qt_ptratio->Graph(part);
    gr_ptratio->SetTitle("N_{#gamma}(p_{T}^{2#gamma})/N_{#gamma}(p_{T}^{1#gamma})");
    aset(gr_ptratio, "p_{T} [GeV]","Ratio", 0.,30., 0.,10.);
    style(gr_ptratio, 20+part, 1+part);
    if(part==0)
      gr_ptratio->Draw("AP");
    else
      gr_ptratio->Draw("P");

    delete h_pt2photon;
    delete h_pt1photon;
  }

  c0->Print("plots/MissingRatio.pdf");
  c1->Print("plots/Merge-photon.pdf");
  c2->Print("plots/PtRatio-sim.pdf");
  qt_miss->Save();
  qt_merge1->Save();
  qt_merge2->Save();
  qt_ptratio->Save();
}
