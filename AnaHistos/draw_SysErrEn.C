#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "SumSysErr.h"

void draw_SysErrEn()
{
  const char *sysname[4] = {"Sum", "Global scale", "Non-lin", "Geom"};

  QueryTree *qt_sys = new QueryTree("data/syserr-en-fast.root", "RECREATE");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-Fast-histo-syserr.root");
  THnSparse *hn_photon = (THnSparse*)f->Get("hn_photon");
  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

  const int s1 = 8;
  const int s2 = 2;
  for(int ipt=4; ipt<npT;ipt+=(ipt<20?s1:s2))
  {
    double ndir[2][4], endir[2][4];  // photon|pion, isys
    double rsys[2][3], ersys[2][3];  // photon|pion, isys-1

    for(int isys=0; isys<4; isys++)
    {
      hn_photon->GetAxis(3)->SetRange(isys+1,isys+1);
      TH1 *h_pt = hn_photon->Projection(1);  // pt_reco
      ndir[0][isys] = h_pt->IntegralAndError(ipt+1,ipt+(ipt<20?s1:s2), endir[0][isys]);
      delete h_pt;

      TH1 *h_minv;
      hn_pion->GetAxis(5)->SetRange(isys+1,isys+1);
      hn_pion->GetAxis(1)->SetRange(ipt+1,ipt+(ipt<20?s1:s2));  // pt_reco
      h_minv = (TH1*)hn_pion->Projection(2);  // minv
      double hpt = ipt<20 ? 0. : 0.01;
      FitMinv(h_minv, ndir[1][isys], endir[1][isys], false, 0.11-hpt,0.16+hpt);
      delete h_minv;

      if(isys)
        for(int itype=0; itype<2; itype++)
        {
          int index = itype + 2*isys;
          double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
          rsys[itype][isys-1] = ndir[itype][isys]/ndir[itype][0];
          ersys[itype][isys-1] = rsys[itype][isys-1]*sqrt(pow(endir[itype][isys]/ndir[itype][isys],2) + pow(endir[itype][0]/ndir[itype][0],2));
          if( TMath::Finite(rsys[itype][isys-1]+ersys[itype][isys-1]) )
            qt_sys->Fill(ipt, index, xpt, fabs(rsys[itype][isys-1]-1), ersys[itype][isys-1]);
        } // itype
    } // isys

    for(int itype=0; itype<2; itype++)
    {
      double sum, esum;
      SumSysErr(3, rsys[itype], ersys[itype], sum, esum);
      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      if( TMath::Finite(sum+esum) )
        qt_sys->Fill(ipt, itype, xpt, sum, esum);
    } // itype
  } // ipt

  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(4);
  qt_sys->Write();
  for(int itype=0; itype<2; itype++)
  {
    mc();
    mcd();
    for(int isys=0; isys<4; isys++)
    {
      int index = itype + 2*isys;
      TGraphErrors *gr = qt_sys->Graph(index);
      aset(gr, "p_{T} [GeV]","SysErr", 5.,30., 0.,0.1);
      style(gr, 20+isys, 1+isys);
      char *opt = isys ? "P" : "AL";
      gr->Draw(opt);
      if(itype == 0)
        leg0->AddEntry(gr, sysname[isys], "P");
    }
    leg0->Draw();
    mcw( 0, Form("sysen-pion%d",itype) );
    char *type = itype ? "pion" : "photon";
    c0->Print(Form("plots/SysErrEn-%s.pdf",type));
  }
  qt_sys->Close();
}

void draw_SysErrEn_Pythia()
{
  const char *sysname[4] = {"Sum", "Global scale", "Non-lin", "Geom"};

  QueryTree *qt_sys = new QueryTree("data/syserr-en-pythia.root", "RECREATE");

  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");
  QueryTree *qt_A = new QueryTree("data/GammaRatio.root");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-PH-histo-syserr.root");
  THnSparse *hn_photon = (THnSparse*)f->Get("hn_photon");
  THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");
  hn_photon->GetAxis(4)->SetRange(3,3); // econe_trk[ival]: EMCal, nomap, withmap
  hn_pion->GetAxis(6)->SetRange(3,3); // econe_trk[ival]: EMCal, nomap, withmap

  for(int ipt=10; ipt<npT; ipt++)
  {
    double ndir[2][4], endir[2][4];  // isolated, isys
    double rsys[2][3], ersys[2][3];  // isolated, isys-1
    double xpt, Miss, eMiss, A, eA;
    qt_miss->Query(ipt, 3, xpt, Miss, eMiss);
    qt_A->Query(ipt, 0, xpt, A, eA);

    for(int isys=0; isys<4; isys++)
    {
      hn_photon->GetAxis(5)->SetRange(isys+1,isys+1);
      hn_pion->GetAxis(7)->SetRange(isys+1,isys+1);

      double nphoton[2], enphoton[2];  // isolated
      for(int isolated=0; isolated<2; isolated++)
      {
        hn_photon->GetAxis(3)->SetRange(isolated+1,2);
        TH1 *h_1photon = hn_photon->Projection(0);  // pt_truth
        nphoton[isolated] = h_1photon->GetBinContent(ipt+1);
        enphoton[isolated] = h_1photon->GetBinError(ipt+1);
        delete h_1photon;
      }

      double npion[2][2], enpion[2][2];  // isoboth, isopair
      for(int isoboth=0; isoboth<2; isoboth++)
        for(int isopair=0; isopair<2; isopair++)
          if(isoboth != 1 || isopair != 1)
          {
            TH1 *h_minv;
            hn_pion->GetAxis(4)->SetRange(isoboth+1,2);  // isoboth
            hn_pion->GetAxis(5)->SetRange(isopair+1,2);  // isopair
            hn_pion->GetAxis(0)->SetRange(ipt+1,ipt+1);  // pt_truth
            h_minv = (TH1*)hn_pion->Projection(2);  // minv
            double hpt = ipt<20 ? 0. : 0.01;
            FitMinv(h_minv, npion[isoboth][isopair], enpion[isoboth][isopair], false, 0.11-hpt,0.16+hpt);
            delete h_minv;
          }

      ndir[0][isys] = nphoton[0] - (1 + Miss)*(1 + A)*npion[0][0];
      endir[0][isys] = sqrt(pow(enphoton[0],2) + pow((1 + Miss)*(1 + A)*enpion[0][0],2));

      ndir[1][isys] = nphoton[1] - npion[1][0] - Miss*(1 + A)*npion[0][1];
      endir[1][isys] = sqrt(pow(enphoton[1],2) + pow(enpion[1][0],2) + pow(Miss*(1 + A)*enpion[0][1],2));

      if(isys)
        for(int isolated=0; isolated<2; isolated++)
        {
          int index = isolated + 2*isys;
          rsys[isolated][isys-1] = ndir[isolated][isys]/ndir[isolated][0];
          ersys[isolated][isys-1] = rsys[isolated][isys-1]*sqrt(pow(endir[isolated][isys]/ndir[isolated][isys],2) + pow(endir[isolated][0]/ndir[isolated][0],2));
          if( TMath::Finite(rsys[isolated][isys-1]+ersys[isolated][isys-1]) )
            qt_sys->Fill(ipt, index, xpt, fabs(rsys[isolated][isys-1]-1), ersys[isolated][isys-1]);
        } // isolated
    } // isys

    for(int isolated=0; isolated<2; isolated++)
    {
      double sum, esum;
      SumSysErr(3, rsys[isolated], ersys[isolated], sum, esum);
      double xpt = (pTbin[ipt] + pTbin[ipt+1]) / 2.;
      if( TMath::Finite(sum+esum) )
        qt_sys->Fill(ipt, isolated, xpt, sum, esum);
    } // isolated
  } // ipt

  const double xmin = 5.;
  const double xmax = 30.;
  TF1 *fn_engl = new TF1("fn_engl", "1+[0]-[1]/([2]+x*x)", xmin,xmax);

  legi(0, 0.2,0.8,0.9,0.9);
  leg0->SetNColumns(4);
  qt_sys->Write();
  for(int isolated=0; isolated<2; isolated++)
  {
    mc();
    mcd();
    for(int isys=0; isys<4; isys++)
    {
      int index = isolated + 2*isys;
      TGraphErrors *gr = qt_sys->Graph(index);
      aset(gr, "p_{T} [GeV]","SysErr", 5.,30., 0.,0.5);
      style(gr, 20+isys, 1+isys);
      char *opt = isys ? "P" : "AL";
      gr->Draw(opt);
      if(isys == 0)
        gr->Fit(fn_engl, "RQ");
      if(isolated == 0)
        leg0->AddEntry(gr, sysname[isys], "P");
    }
    leg0->Draw();
    mcw( 0, Form("sysen-iso%d",isolated) );
    char *type = isolated ? "iso" : "inc";
    c0->Print(Form("plots/SysErrEn-%sphoton.pdf",type));
  }
  qt_sys->Close();
}
