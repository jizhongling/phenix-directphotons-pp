#include "DivideFunctions.h"

void draw_Smear_Sasha()
{
  TFile *f1 = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo.root");
  THnSparse *hn_pion = (THnSparse*)f1->Get("hn_pion");

  TFile *f2 = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/pi0cross_run13pp510gev/fastMC/eff-histo.root");
  THnSparse *hn_pion2 = (THnSparse*)f2->Get("hn_pion");

  TFile *f[3];
  TH1 *hpt_acc1[3];
  TH1 *hpt_acc2[3];

  mc(0);
  mc(1);

  legi(0, 0.2,0.8,0.9,0.9);
  legi(1, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);
  leg1->SetNColumns(3);

  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};
  const char *pname[3] = {"PbScW", "PbScE", "PbGlE"};  // must be non-const if used directly in TLegend

  for(int part=0; part<3; part++)
  {
    hn_pion->GetAxis(3)->SetRange(secl[part],sech[part]);
    TH1 *h_pt0 = hn_pion->Projection(0);
    TH1 *h_pt1 = hn_pion->Projection(1);
    TGraphErrors *gr1 = DivideHisto(h_pt1, h_pt0);
    delete h_pt0;
    delete h_pt1;

    //hn_pion2->GetAxis(3)->SetRange(secl[part],sech[part]);
    //TH1 *h_pt0 = hn_pion2->Projection(0);
    //TH1 *h_pt1 = hn_pion2->Projection(1);
    //TGraphErrors *gr2 = DivideHisto(h_pt1, h_pt0);
    //delete h_pt0;
    //delete h_pt1;

    f[part] = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/pi0cross_run13pp510gev/fastMC/eff-%d.root",part));
    hpt_acc1[part] = (TH1*)f[part]->Get("hpt_acc1");
    hpt_acc2[part] = (TH1*)f[part]->Get("hpt_acc2");
    TGraphErrors *gr2 = DivideHisto(hpt_acc2[part], hpt_acc1[part]);
    //TGraphAsymmErrors *gr2 = new TGraphAsymmErrors(hpt_acc2[part], hpt_acc1[part], "n");

    mcd(0);
    TGraphErrors *gr = DivideGraph(gr1, gr2);
    gr->SetTitle("Smear Ratio");
    aset(gr, "p_{T} (GeV/c)", "#frac{My smear}{Sasha's smear}", 0.,30., 0.7,1.3);
    style(gr, 20+part, 1+part);
    if(part==0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    leg0->AddEntry(gr, pname[part], "P");
    leg0->Draw();

    mcd(1);
    gr1->SetTitle("Smear");
    aset(gr1, "p_{T} (GeV/c)", "Smear", 0.,30., 0.7,1.5);
    style(gr1, 20+part, 1+part);
    style(gr2, 20+part, 1+part);
    gr2->SetLineWidth(6.);
    if(part==0)
    {
      gr1->Draw("AP");
      gr2->Draw("L");
    }
    else
    {
      gr1->Draw("P");
      gr2->Draw("L");
    }
    leg1->AddEntry(gr1, pname[part], "P");
    leg1->Draw();

    gr2->RemovePoint(0);
    for(int ip=1; ip<3; ip++)
    {
      gr1->RemovePoint(ip);
      gr2->RemovePoint(ip);
      gr->RemovePoint(ip);
    }
  }

  c0->Print("plots/SmearRatio.pdf");
  c1->Print("plots/Smear.pdf");
}
