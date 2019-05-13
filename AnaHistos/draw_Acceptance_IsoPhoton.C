#include "GlobalVars.h"
#include "QueryTree.h"
#include "DivideFunctions.h"

void draw_Acceptance_IsoPhoton()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  QueryTree *qt_acc = new QueryTree("data/Acceptance-isophoton.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo.root");
  TH1 *h_photon = (TH1*)f->Get("h_photon");
  TH1 *h_isophoton = (TH1*)f->Get("h_isophoton");
  TH1 *h_isolated = (TH1*)f->Get("h_isolated");
  THnSparse *hn_geom = (THnSparse*)f->Get("hn_geom");
  THnSparse *hn_isolated = (THnSparse*)f->Get("hn_isolated");

  TGraphErrors *gr_geom[3];
  TGraphErrors *gr_iso[3];
  TGraphErrors *gr_smear[3];

  for(int part=0; part<3; part++)
  {
    hn_geom->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_isolated->GetAxis(2)->SetRange(secl[part],sech[part]);
    TH1 *h_geom = hn_geom->Projection(0);
    TH1 *h_iso = hn_isolated->Projection(0);
    TH1 *h_reco = hn_isolated->Projection(1);

    gr_geom[part] = DivideHisto(h_geom, h_photon);
    gr_iso[part] = DivideHisto(h_iso, h_geom);
    gr_smear[part] = DivideHisto(h_reco, h_iso);
    qt_acc->Fill(h_reco, h_isolated, part);

    delete h_geom;
    delete h_iso;
    delete h_reco;
  }

  for(int ic=0; ic<4; ic++)
    mc(ic);
  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(3);

  for(int part=0; part<3; part++)
  {
    mcd(0);
    gr_geom[part]->SetTitle("Geometric acceptance");
    aset(gr_geom[part], "p_{T} [GeV]","InAcc/All", 4.,30., 0.,0.12);
    style(gr_geom[part], part+20, part+1);
    if(part==0)
      gr_geom[part]->Draw("AP");
    else
      gr_geom[part]->Draw("P");

    mcd(1);
    gr_iso[part]->SetTitle("Isolated over inclusive prompt photons");
    aset(gr_iso[part], "p_{T} [GeV]","Iso/InAcc", 4.,30.);
    style(gr_iso[part], part+20, part+1);
    if(part==0)
      gr_iso[part]->Draw("AP");
    else
      gr_iso[part]->Draw("P");

    mcd(2);
    gr_smear[part]->SetTitle("Photon p_{T} smearing");
    aset(gr_smear[part], "p_{T} [GeV]","Reco/Truth", 4.,30., 0.9,1.5);
    style(gr_smear[part], part+20, part+1);
    if(part==0)
      gr_smear[part]->Draw("AP");
    else
      gr_smear[part]->Draw("P");

    mcd(3);
    TGraphErrors *gr_acc = qt_acc->Graph(part);
    gr_acc->SetTitle("Combined acceptance");
    aset(gr_acc, "p_{T} [GeV]","Acceptance", 4.,30., 0.,0.4);
    style(gr_acc, part+20, part+1);
    if(part==0)
      gr_acc->Draw("AP");
    else
      gr_acc->Draw("P");
    leg0->AddEntry(gr_acc, pname[part], "P");
  }

  leg0->Draw();
  qt_acc->Save();
  c3->Print("plots/Acceptance-isophoton.pdf");
}
