#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_CrossSection_IsoPhoton()
{
  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  const double PI = TMath::Pi();
  const double DeltaEta = 0.5;
  //const double NBBC =  3.48e11;  // from DAQ
  const double NBBC =  3.52e11;  // from rejection power
  const double XBBC = 32.51e9;
  const double eXBBC = 3.24e9;
  const double Pile[3] = {1.10, 1.08, 1.05};
  const double ePile = 0.02;
  const double TrigBBC = 0.91;
  const double eTrigBBC = 0.01;
  const double Prob = 0.96;
  const double eProb = 0.02;
  const double ToF[3] = {0.992, 0.992, 0.997};
  const double eToF[3] = {0.002, 0.002, 0.002};
  const double Conv[3] = {0.849, 0.971, 0.959};
  const double eConv[3] = {0.027, 0.023, 0.029};
  const double Norm[3] = {0.320, 0.321, 0.250};
  const double eNorm[3] = {0.005, 0.007, 0.005};
  const double A = 0.28;
  const double eA = 0.05;

  // function for pT weight for direct photon
  cross_ph = new TF1("cross_ph", "x**(-[1]-[2]*log(x/[0]))*(1-(x/[0])**2)**[3]*[4]", 0, 30);
  cross_ph->SetParameters(255., 5.98, 0.273, 14.43, 1.);
  double ndata = 0.;
  double nfit = 0.;

  QueryTree *qt_cross = new QueryTree("data/CrossSection-isophoton.root", "RECREATE");

  QueryTree *qt_acc = new QueryTree("data/Acceptance-isophoton.root");
  QueryTree *qt_ert = new QueryTree("data/ERTEff-photon.root");
  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");
  QueryTree *qt_miss_eta = new QueryTree("data/MissingRatio-eta.root");
  QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root");
  QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root");
  QueryTree *qt_badpass = new QueryTree("data/MergePassRate.root");
  QueryTree *qt_veto = new QueryTree("data/SelfVeto.root");
  QueryTree *qt_sys = new QueryTree("data/syserr-en-fast.root");
  QueryTree *qt_rbg = new QueryTree("data/BgRatio.root");
  QueryTree *qt_pt = new QueryTree("data/PtShift.root");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[evtype][part]
  TH1 *h_1photon[3][3];
  TH2 *h2_isoboth[3][3];
  TH2 *h2_isopair[3][3];
  TH2 *h2_isopair2pt[3][3];

  int checkmap = 1;
  int ival = 1;

  TH1 *h_1photon_t = (TH1*)f->Get("h_1photon_0");
  h_1photon_t = (TH1*)h_1photon_t->Clone();
  h_1photon_t->Reset();

  TH2 *h2_2photon_t = (TH2*)f->Get("h2_2photon_0");
  h2_2photon_t = (TH2*)h2_2photon_t->Clone();
  h2_2photon_t->Reset();

  for(int evtype=0; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h_1photon[evtype][part] = (TH1*)h_1photon_t->Clone(Form("h_1photon_type%d_part%d",evtype,part));
      h2_isoboth[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isoboth_type%d_part%d",evtype,part));
      h2_isopair[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isopair_type%d_part%d",evtype,part));
      h2_isopair2pt[evtype][part] = (TH2*)h2_2photon_t->Clone(Form("h2_isopair2pt_type%d_part%d",evtype,part));

      int isolated = 1;
      int ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isolated + 3*3*2*2*ival;
      TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_%d",ih));
      h_1photon[evtype][part]->Add(h_tmp);

      for(int isopair=0; isopair<2; isopair++)
      {
        int isoboth = 1;
        int ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isoboth + 3*3*2*2*isopair + 3*3*2*2*2*ival;
        TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
        h2_isoboth[evtype][part]->Add(h2_tmp);
      } // isopair

      for(int isoboth=0; isoboth<2; isoboth++)
      {
        int isopair = 1;
        int ih = part + 3*evtype + 3*3*checkmap + 3*3*2*isoboth + 3*3*2*2*isopair + 3*3*2*2*2*ival;
        TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon_%d",ih));
        h2_isopair[evtype][part]->Add(h2_tmp);

        h2_tmp = (TH2*)f->Get(Form("h2_2photon2pt_%d",ih));
        h2_isopair2pt[evtype][part]->Add(h2_tmp);
      } // isoboth
    }

  for(int part=0; part<3; part++)
  {
    mc(part, 6,5);
    mc(part+3, 6,5);
    mc(part+6, 6,5);
  }

  for(int ipt=2; ipt<npT; ipt++)
  {
    double xpt, yy[3][3], eyy[3][3];  // isys, part

    double SysPhoton, eSysPhoton, SysPion, eSysPion;
    qt_sys->Query(ipt<20?ipt/4*4:ipt, 0, xpt, SysPhoton, eSysPhoton);
    qt_sys->Query(ipt<20?ipt/4*4:ipt, 1, xpt, SysPion, eSysPion);

    for(int part=0; part<3; part++)
    {
      int evtype = 2;
      if(ipt < 22)  // <14GeV use ERT_4x4c
        evtype = 2;
      else  // >14GeV use ERT_4x4b
        evtype = 1;

      double xpt, Acc, eAcc, TrigERT, eTrigERT, Miss, eMiss, MissEta, eMissEta,
             Merge1, eMerge1, Merge2, eMerge2, BadPass, eBadPass, Veto, eVeto, rbg, erbg;
      qt_rbg->Query(ipt, part/2, xpt, rbg, erbg);
      qt_acc->Query(ipt, part, xpt, Acc, eAcc);
      qt_ert->Query(ipt, part, xpt, TrigERT, eTrigERT);
      qt_miss->Query(ipt, part, xpt, Miss, eMiss);
      qt_miss_eta->Query(ipt, part, xpt, MissEta, eMissEta);
      qt_veto->Query(ipt, part, xpt, Veto, eVeto);
      if(xpt > 12.)
      {
        qt_merge1->Query(ipt, part, xpt, Merge1, eMerge1);
        qt_merge2->Query(ipt, part, xpt, Merge2, eMerge2);
        qt_badpass->Query(ipt, part/2, xpt, BadPass, eBadPass);
      }
      else
      {
        Merge1 = eMerge1 = Merge2 = eMerge2 = BadPass = eBadPass = 0.;
      }
      if(!TMath::Finite(rbg+erbg) || xpt > 16.)
      {
        rbg = erbg = 0.1;
      }

      if(ipt >= 20)
      {
        if(part == 0)
        {
          TrigERT = 0.971;
          eTrigERT = 0.003;
        }
        else if(part == 1)
        {
          TrigERT = 0.950;
          eTrigERT = 0.004;
        }
        else
        {
          TrigERT = 0.693;
          eTrigERT = 0.006;
        }
      }

      double nphoton = h_1photon[evtype][part]->GetBinContent(ipt+1);
      double enphoton = sqrt(nphoton);

      TH1 *h_minv;

      mcd(part, ipt+1);
      double nisoboth = 1., enisoboth = 1.;
      h_minv = (TH1*)h2_isoboth[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, nisoboth, enisoboth, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, nisoboth, enisoboth, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, nisoboth, enisoboth, false, 0.10,0.17);
      nisoboth /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      mcd(part+3, ipt+1);
      double nisopair = 1., enisopair = 1.;
      h_minv = (TH1*)h2_isopair[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, nisopair, enisopair, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, nisopair, enisopair, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, nisopair, enisopair, false, 0.10,0.17);
      nisopair /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      mcd(part+6, ipt+1);
      double nisopair2pt = 1., enisopair2pt = 1.;
      h_minv = (TH1*)h2_isopair2pt[evtype][part]->ProjectionY("h_py", ipt+1,ipt+1)->Clone("h_minv");
      h_minv->Rebin(10);
      h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin[ipt], pTbin[ipt+1]) );
      if(ipt < 20)  // <10GeV +-25MeV; >10GeV +-35MeV
        FitMinv(h_minv, nisopair2pt, enisopair2pt, true, 0.11,0.16);
      else if(ipt < 23)  // <16GeV subtract background
        FitMinv(h_minv, nisopair2pt, enisopair2pt, true, 0.10,0.17);
      else  // >16GeV don't subtract background
        FitMinv(h_minv, nisopair2pt, enisopair2pt, false, 0.10,0.17);
      nisopair2pt /= bck[part/2][ipt] * meff[part/2][ipt];
      delete h_minv;

      double Eff = Conv[part] * Prob * ToF[part];
      double AIso = A * (Veto+MissEta)/(1+2.*MissEta) * (1+2.*Miss+Merge1);
      double nbg = (1 + Merge1*Conv[part]*(1-Conv[part])) * nisoboth/Eff/Eff + Miss * nisopair/Eff/Eff + Merge2/2.*BadPass * nisopair2pt/Eff/Eff + AIso * nisopair/Eff/Eff;

      for(int isys=0; isys<3; isys++)
      {
        double ndir, erel;
        if(isys < 2)
        {
          ndir = nphoton/Eff - nbg*(1 + rbg*isys);
          erel = 1e-9;
        }
        else
        {
          ndir = nphoton/Eff*(1 + SysPhoton) - nbg*(1 + (ipt>8?1:-1)*SysPion);
          erel = sqrt(pow(nphoton/Eff*eSysPhoton,2) + pow(nbg*eSysPion,2));
          erel /= ndir;
        }
        double endir = sqrt((pow(enisoboth,2)*pow(1 + (1 - Conv[part])*Conv[part]*Merge1,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (0.25*pow(BadPass,2)*pow(enisopair2pt,2)*pow(Merge2,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (pow(A,2)*pow(eVeto,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(nisopair,2))/(pow(Conv[part],4)*pow(1 + 2.*MissEta,2)*pow(Prob,4)*pow(ToF[part],4)) + (0.25*pow(BadPass,2)*pow(eMerge2,2)*pow(nisopair2pt,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + (0.25*pow(eBadPass,2)*pow(Merge2,2)*pow(nisopair2pt,2))/(pow(Conv[part],4)*pow(Prob,4)*pow(ToF[part],4)) + pow(enphoton,2)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2)) + (pow(eA,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(nisopair,2)*pow(MissEta + Veto,2))/(pow(Conv[part],4)*pow(1 + 2.*MissEta,2)*pow(Prob,4)*pow(ToF[part],4)) + pow(eToF[part],2)*pow((2*(1 + (1 - Conv[part])*Conv[part]*Merge1)*nisoboth)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],3)) + (2*Miss*nisopair)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],3)) + (1.*BadPass*Merge2*nisopair2pt)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],3)) - nphoton/(Conv[part]*Prob*pow(ToF[part],2)) + (2*A*(1 + Merge1 + 2.*Miss)*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],3)),2) + pow(eProb,2)*pow((2*(1 + (1 - Conv[part])*Conv[part]*Merge1)*nisoboth)/(pow(Conv[part],2)*pow(Prob,3)*pow(ToF[part],2)) + (2*Miss*nisopair)/(pow(Conv[part],2)*pow(Prob,3)*pow(ToF[part],2)) + (1.*BadPass*Merge2*nisopair2pt)/(pow(Conv[part],2)*pow(Prob,3)*pow(ToF[part],2)) - nphoton/(Conv[part]*pow(Prob,2)*ToF[part]) + (2*A*(1 + Merge1 + 2.*Miss)*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,3)*pow(ToF[part],2)),2) + pow(enisopair,2)*pow(-(Miss/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2))) - (A*(1 + Merge1 + 2.*Miss)*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2)),2) + pow(eMissEta,2)*pow(-((A*(1 + Merge1 + 2.*Miss)*nisopair)/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2))) + (2.*A*(1 + Merge1 + 2.*Miss)*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*pow(1 + 2.*MissEta,2)*pow(Prob,2)*pow(ToF[part],2)),2) + pow(eMiss,2)*pow(-(nisopair/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2))) - (2.*A*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2)),2) + pow(eMerge1,2)*pow(-(((1 - Conv[part])*nisoboth)/(Conv[part]*pow(Prob,2)*pow(ToF[part],2))) - (A*nisopair*(MissEta + Veto))/(pow(Conv[part],2)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2)),2) + pow(eConv[part],2)*pow(-((((1 - Conv[part])*Merge1 - Conv[part]*Merge1)*nisoboth)/(pow(Conv[part],2)*pow(Prob,2)*pow(ToF[part],2))) + (2*(1 + (1 - Conv[part])*Conv[part]*Merge1)*nisoboth)/(pow(Conv[part],3)*pow(Prob,2)*pow(ToF[part],2)) + (2*Miss*nisopair)/(pow(Conv[part],3)*pow(Prob,2)*pow(ToF[part],2)) + (1.*BadPass*Merge2*nisopair2pt)/(pow(Conv[part],3)*pow(Prob,2)*pow(ToF[part],2)) - nphoton/(pow(Conv[part],2)*Prob*ToF[part]) + (2*A*(1 + Merge1 + 2.*Miss)*nisopair*(MissEta + Veto))/(pow(Conv[part],3)*(1 + 2.*MissEta)*pow(Prob,2)*pow(ToF[part],2)),2));

        if(ipt >= 22)  // >14GeV use ERT_4x4b
        {
          ndir *= Norm[part];
          endir *= Norm[part];
        }

        double dummy;
        qt_pt->Query(ipt, 0, dummy, xpt, dummy);
        yy[isys][part] = (XBBC/NBBC) / (2*PI*xpt) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
          * ndir / Acc / TrigERT / TrigBBC * Pile[part];
        eyy[isys][part] = yy[isys][part] * sqrt( pow(endir/ndir,2)
            + pow(eAcc/Acc,2)
            + pow(eTrigERT/TrigERT,2)
            + pow(ePile/Pile[part],2)
            //+ pow(eTrigBBC/TrigBBC,2) + pow(eXBBC/XBBC,2)
            );
        if(isys)
        {
          double rsys = yy[isys][part]/yy[0][part];
          double ersys = rsys*erel;
          if( TMath::Finite(rsys+ersys) )
            qt_cross->Fill(ipt, 1+part+3*isys, xpt, fabs(rsys-1), ersys);
        }
      } // isys
      if( TMath::Finite(yy[0][part]+eyy[0][part]) )
        qt_cross->Fill(ipt, part, xpt, yy[0][part], eyy[0][part]);
    } // part

    double ybar, eybar;
    Chi2Fit(3, yy[0], eyy[0], ybar, eybar);
    if( TMath::Finite(ybar + eybar) )
      qt_cross->Fill(ipt, 3, xpt, ybar, eybar);
    if( ipt >= 12 )
    {
      ndata += ybar;
      nfit += cross_ph->Eval(xpt);
    }
  } // ipt
  cross_ph->SetParameter(4, ndata/nfit);

  mc(9, 2,1);
  legi(0, 0.4,0.7,0.7,0.9);

  for(int part=0; part<4; part++)
  {
    TGraphErrors *gr = qt_cross->Graph(part);
    mcd(9, part/3+1);
    gPad->SetLogy();
    if(part == 0)
      gr->SetTitle("Separated");
    else if(part == 3)
      gr->SetTitle("Combined");
    aset(gr, "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{-3}]", 6.1,30., 1e-1, 2e3);
    style(gr, part+20, part+1);
    if(part%3==0)
    {
      gr->Draw("AP");
      cross_ph->DrawCopy("LSAME");
    }
    else
      gr->Draw("P");
    if(part<3)
      leg0->AddEntry(gr, pname[part], "P");
  }

  mcd(9, 1);
  leg0->Draw();
  mcd(9, 2);
  c9->Print("plots/CrossSection-isophoton.pdf");

  mc(10);
  mcd(10);
  legi(1, 0.2,0.8,0.9,0.9);
  leg1->SetNColumns(3);

  for(int isys=1; isys<3; isys++)
  {
    for(int part=0; part<3; part++)
    {
      TGraphErrors *gr = qt_cross->Graph(1+part+3*isys);
      aset(gr, "p_{T} [GeV]", "SysErr", 6.1,30., 0.,0.06+0.14*(isys-1));
      style(gr, part+20, part+1);
      char *opt = part==0 ? "AP" : "P";
      gr->Draw(opt);
      if(isys == 1)
        leg1->AddEntry(gr, pname[part], "P");
    }
    leg1->Draw();
    char *type = isys==1 ? "Fit" : "En";
    c10->Print(Form("plots/SysErr%s-isophoton.pdf",type));
  }

  qt_cross->Write();
  for(int part=0; part<3; part++)
  {
    mcw( part, Form("Minv-isoboth-part%d",part) );
    mcw( part+3, Form("Minv-isopair-part%d",part) );
    mcw( part+6, Form("Minv-isopair2pt-part%d",part) );
  }
  qt_cross->Close();
}
