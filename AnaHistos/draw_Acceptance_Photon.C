#include "GlobalVars.h"
#include "DivideFunctions.h"

void draw_Acceptance_Photon()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  TGraphAsymmErrors *gr_geom[3];
  TGraphAsymmErrors *gr_iso[3];
  TGraphErrors *gr_smear[3];
  TGraphAsymmErrors *gr_acc[3];
  for(int part=0; part<3; part++)
  {
    gr_geom[part] = new TGraphAsymmErrors(npT);
    gr_iso[part] = new TGraphAsymmErrors(npT);
    gr_acc[part] = new TGraphAsymmErrors(npT);
    gr_acc[part]->SetName(Form("gr_%d",part));
  }

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo.root");
  TH1 *h_photon = (TH1*)f->Get("h_photon");
  TH1 *h_isophoton = (TH1*)f->Get("h_isophoton");
  TH1 *h_isolated = (TH1*)f->Get("h_isolated");
  THnSparse *hn_geom = (THnSparse*)f->Get("hn_geom");
  THnSparse *hn_isolated = (THnSparse*)f->Get("hn_isolated");

  for(int part=0; part<3; part++)
  {
    hn_geom->GetAxis(2)->SetRange(secl[part],sech[part]);
    hn_isolated->GetAxis(2)->SetRange(secl[part],sech[part]);
    TH1 *h_geom = hn_geom->Projection(0);
    TH1 *h_iso = hn_isolated->Projection(0);
    TH1 *h_reco = hn_isolated->Projection(1);
    for(int binx=1; binx<=h_photon->GetNbinsX(); binx++)
      if( h_reco->GetBinContent(binx) > h_isolated->GetBinContent(binx) )
        h_reco->SetBinContent(binx, h_isolated->GetBinContent(binx));

    gr_geom[part]->Divide(h_geom, h_photon, "n");
    gr_iso[part]->Divide(h_iso, h_geom, "n");
    gr_smear[part] = DivideHisto(h_reco, h_iso);
    gr_acc[part]->Divide(h_reco, h_isolated, "n");
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
    aset(gr_geom[part], "p_{T} [GeV]","Incl/All", 4.,30., 0.,0.12);
    style(gr_geom[part], part+20, part+1);
    if(part==0)
      gr_geom[part]->Draw("AP");
    else
      gr_geom[part]->Draw("P");

    mcd(1);
    gr_iso[part]->SetTitle("Isolated over inclusive photons");
    aset(gr_iso[part], "p_{T} [GeV]","Iso/Incl", 4.,30.);
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
    gr_acc[part]->SetTitle("Combined acceptance");
    aset(gr_acc[part], "p_{T} [GeV]","Acceptance", 4.,30., 0.,0.4);
    style(gr_acc[part], part+20, part+1);
    if(part==0)
      gr_acc[part]->Draw("AP");
    else
      gr_acc[part]->Draw("P");
    leg0->AddEntry(gr_acc[part], pname[part], "P");
  }

  leg0->Draw();

  TFile *f_out = new TFile("data/Acceptance-photon.root", "RECREATE");
  for(int part=0; part<3; part++)
  {
    gr_acc[part]->Write();
  }
  f_out->Close();
}
