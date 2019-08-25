#include "GlobalVars.h"
#include "QueryTree.h"
#include "DivideFunctions.h"
#include "FitMinv.h"

void draw_Acceptance_IsoPhoton()
{
  const char *simname[2] = {"FastMC", "PISA"};
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  const double Conv[3] = {0.8855, 0.9913, 0.9913};
  const double eConv[3] = {2e-4, 7e-5, 7e-5};
  const double A = 0.24;
  const double eA = 0.04;

  const double Prob[2][npT] = {
    { 1, 1, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 },
    { 1, 1, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97 }
  };
  const double eProb = 0.02;

  QueryTree *qt_acc = new QueryTree("data/Acceptance-isophoton.root", "RECREATE");

  //QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");
  //QueryTree *qt_miss_eta = new QueryTree("data/MissingRatio-eta.root");
  //QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root");
  //QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root");
  //QueryTree *qt_badpass = new QueryTree("data/MergePassRate.root");
  //QueryTree *qt_veto = new QueryTree("data/SelfVeto.root");

  TFile *f_pythia = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-GenPH-histo.root");
  TH1 *h_photon = (TH1*)f_pythia->Get("h_photon_eta025");
  TH1 *h_isolated = (TH1*)f_pythia->Get("h_isophoton_eta025");
  THnSparse *hn_geom = (THnSparse*)f_pythia->Get("hn_geom");
  THnSparse *hn_isolated = (THnSparse*)f_pythia->Get("hn_isolated");

  TFile *f_pisa = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/HadronResponse-histo.root");
  THnSparse *hn_1photon = (THnSparse*)f_pisa->Get("hn_1photon");
  THnSparse *hn_2photon = (THnSparse*)f_pisa->Get("hn_2photon");

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

  for(int ic=0; ic<3; ic++)
    mc(ic);
  mc(3, 2,1);
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

    //mc(4+part, 6,5);
    //mc(4+part+3, 6,5);
    //mc(4+part+6, 6,5);

    for(int ipt=10; ipt<npT; ipt++)
    {
      double nisophoton = h_isolated->GetBinContent(ipt+1);
      double enisophoton = h_isolated->GetBinError(ipt+1);

      hn_1photon->GetAxis(2)->SetRange(1,2);  // isolated
      hn_1photon->GetAxis(1)->SetRange(secl[part],sech[part]);  // sector
      TH1 *h_1photon = hn_1photon->Projection(0);  // pt
      double nphoton = h_1photon->GetBinContent(ipt+1);
      double enphoton = h_1photon->GetBinError(ipt+1);
      delete h_1photon;

      //TH1 *h_minv;

      //mcd(4+part, ipt+1);
      //double nisoboth = 1., enisoboth = 1.;
      //hn_2photon->GetAxis(4)->SetRange(2,2);  // isoboth
      //hn_2photon->GetAxis(5)->SetRange(1,2);  // isopair
      //hn_2photon->GetAxis(3)->SetRange(secl[part],sech[part]);  // sector
      //hn_2photon->GetAxis(0)->SetRange(ipt+1,ipt+1);  // pt_1photon
      //hn_2photon->GetAxis(1)->SetRange(1,npT);  // pt_2photon
      //h_minv = (TH1*)hn_2photon->Projection(2);  // minv
      //h_minv->Rebin(10);
      //h_minv->SetTitle( Form("part %d p^{isoboth}_{T}: %3.1f-%3.1f GeV", part, pTbin[ipt], pTbin[ipt+1]) );
      //// don't subtract background
      //FitMinv(h_minv, nisoboth, enisoboth, false, 0.10,0.17);
      //nisoboth /= 1.1;
      //if( !TMath::Finite(enisoboth) )
      //  enisoboth = sqrt(nisoboth);
      //delete h_minv;

      //mcd(4+part+3, ipt+1);
      //double nisopair = 1., enisopair = 1.;
      //hn_2photon->GetAxis(4)->SetRange(1,2);  // isoboth
      //hn_2photon->GetAxis(5)->SetRange(2,2);  // isopair
      //hn_2photon->GetAxis(3)->SetRange(secl[part],sech[part]);  // sector
      //hn_2photon->GetAxis(0)->SetRange(ipt+1,ipt+1);  // pt_1photon
      //hn_2photon->GetAxis(1)->SetRange(1,npT);  // pt_2photon
      //h_minv = (TH1*)hn_2photon->Projection(2);  // minv
      //h_minv->Rebin(10);
      //h_minv->SetTitle( Form("part %d p^{isopair}_{T}: %3.1f-%3.1f GeV", part, pTbin[ipt], pTbin[ipt+1]) );
      //if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
      //  FitMinv(h_minv, nisopair, enisopair, true, 0.11,0.16);
      //else if(ipt < 23)  // <16GeV subtract background
      //  FitMinv(h_minv, nisopair, enisopair, true, 0.10,0.17);
      //else  // >16GeV don't subtract background
      //  FitMinv(h_minv, nisopair, enisopair, false, 0.10,0.17);
      //nisopair /= bck[part/2][ipt] * meff[part/2][ipt];
      //if( !TMath::Finite(enisopair) )
      //  enisopair = sqrt(nisopair);
      //delete h_minv;

      //mcd(4+part+6, ipt+1);
      //double nisopair2pt = 1., enisopair2pt = 1.;
      //hn_2photon->GetAxis(4)->SetRange(1,2);  // isoboth
      //hn_2photon->GetAxis(5)->SetRange(2,2);  // isopair
      //hn_2photon->GetAxis(3)->SetRange(secl[part],sech[part]);  // sector
      //hn_2photon->GetAxis(0)->SetRange(1,npT);  // pt_1photon
      //hn_2photon->GetAxis(1)->SetRange(ipt+1,ipt+1);  // pt_2photon
      //h_minv = (TH1*)hn_2photon->Projection(2);  // minv
      //h_minv->Rebin(10);
      //h_minv->SetTitle( Form("part %d p^{isopair2pt}_{T}: %3.1f-%3.1f GeV", part, pTbin[ipt], pTbin[ipt+1]) );
      //if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
      //  FitMinv(h_minv, nisopair2pt, enisopair2pt, true, 0.11,0.16);
      //else if(ipt < 23)  // <16GeV subtract background
      //  FitMinv(h_minv, nisopair2pt, enisopair2pt, true, 0.10,0.17);
      //else  // >16GeV don't subtract background
      //  FitMinv(h_minv, nisopair2pt, enisopair2pt, false, 0.10,0.17);
      //nisopair2pt /= bck[part/2][ipt] * meff[part/2][ipt];
      //if( !TMath::Finite(enisopair2pt) )
      //  enisopair2pt = sqrt(nisopair2pt);
      //delete h_minv;

      //double xpt, Miss, eMiss, MissEta, eMissEta, Merge1, eMerge1, Merge2, eMerge2, BadPass, eBadPass, Veto, eVeto;
      //qt_miss->Query(ipt, part, xpt, Miss, eMiss);
      //qt_miss_eta->Query(ipt, part, xpt, MissEta, eMissEta);
      //qt_merge1->Query(ipt, part, xpt, Merge1, eMerge1);
      //qt_merge2->Query(ipt, part, xpt, Merge2, eMerge2);
      //qt_badpass->Query(ipt, part/2, xpt, BadPass, eBadPass);
      //qt_veto->Query(ipt, part, xpt, Veto, eVeto);

      //double AIso = A * Veto * (1.+MissEta)/(1.+2.*MissEta) * (1+2.*Miss+Merge1);
      //double ndir = nphoton/Conv[part] - (1. + Merge1*Conv[part]*(1.-Conv[part])) * nisoboth/pow(Conv[part],2) - Miss * nisopair/pow(Conv[part],2) - Merge2/2.*BadPass * nisopair2pt - AIso * nisopair;
      //double endir = sqrt(pow(enphoton,2)/pow(Conv[part],2) + (pow(enisoboth,2)* pow(1. + (1. - Conv[part])*Conv[part]*Merge1,2))/ pow(Conv[part],4) + 0.25*pow(BadPass,2)*pow(enisopair2pt,2)* pow(Merge2,2) + 0.25*pow(BadPass,2)*pow(eMerge2,2)* pow(nisopair2pt,2) + 0.25*pow(eBadPass,2)*pow(Merge2,2)* pow(nisopair2pt,2) + (pow(A,2)*pow(eVeto,2)* pow(1. + Miss,2)* pow(1 + Merge1 + 2.*Miss,2)* pow(nisopair,2))/pow(1. + 2.*Miss,2) + pow(eConv[part],2)* pow(-((((1. - Conv[part])*Merge1 - Conv[part]*Merge1)* nisoboth)/pow(Conv[part],2)) + (2*(1. + (1. - Conv[part])*Conv[part]*Merge1)* nisoboth)/pow(Conv[part],3) + (2*Miss*nisopair)/pow(Conv[part],3) - nphoton/pow(Conv[part],2),2) + (pow(eA,2)*pow(1. + Miss,2)* pow(1 + Merge1 + 2.*Miss,2)* pow(nisopair,2)*pow(Veto,2))/ pow(1. + 2.*Miss,2) + pow(enisopair,2)* pow(-(Miss/pow(Conv[part],2)) - (A*(1. + Miss)*(1 + Merge1 + 2.*Miss)* Veto)/(1. + 2.*Miss),2) + pow(eMerge1,2)* pow(-(((1. - Conv[part])*nisoboth)/Conv[part]) - (A*(1. + Miss)*nisopair*Veto)/ (1. + 2.*Miss),2) + pow(eMiss,2)*pow(-(nisopair/ pow(Conv[part],2)) - (2.*A*(1. + Miss)*nisopair*Veto)/ (1. + 2.*Miss) + (2.*A*(1. + Miss)*(1 + Merge1 + 2.*Miss)* nisopair*Veto)/pow(1. + 2.*Miss,2) - (A*(1 + Merge1 + 2.*Miss)*nisopair*Veto)/ (1. + 2.*Miss),2));

      double xpt = ( pTbin[ipt] + pTbin[ipt+1] ) / 2.;
      double ndir = nphoton / Prob[part/2][ipt] / Conv[part];
      double endir = ndir * sqrt( pow(enphoton/nphoton,2) + pow(eProb/Prob[part/2][ipt],2) + pow(eConv[part]/Conv[part],2) );
      double Acc = ndir / nisophoton;
      double eAcc = Acc * sqrt( pow(enisophoton/nisophoton,2) + pow(endir/ndir,2) );
      if( TMath::Finite(Acc+eAcc) )
        qt_acc->Fill(ipt, part+3, xpt, Acc, eAcc);
      else
        cout << "Wrong!!! Not Finite!!!" << endl;
    } // ipt

    for(int pisa=0; pisa<2; pisa++)
    {
      mcd(3, pisa+1);
      TGraphErrors *gr_acc = qt_acc->Graph(part + 3*pisa);
      gr_acc->SetTitle( Form("Combined acceptance in %s",simname[pisa]) );
      aset(gr_acc, "p_{T} [GeV]","Acceptance", 4.,30., 0.,0.4);
      style(gr_acc, part+24, part+1);
      if(part==0)
        gr_acc->Draw("AP");
      else
        gr_acc->Draw("P");
      if(pisa==1)
        leg0->AddEntry(gr_acc, pname[part], "P");
    } // pisa
  } // part

  leg0->Draw();
  c3->Print("plots/Acceptance-isophoton.pdf");
  qt_acc->cd();
  for(int ic=0; ic<4; ic++)
    mcw( ic, Form("c%d",ic) );
  qt_acc->Write();
}
