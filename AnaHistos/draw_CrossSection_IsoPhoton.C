#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"
#include "Chi2Fit.h"

void draw_CrossSection_IsoPhoton()
{
  const int nsys = 18;
  const int ngsys = 5;
  const int gsys[ngsys+1] = {1, 3, 4, 5, 7, nsys+1};
  const char *sysname[ngsys+1] = {"Sum", "En", "Fit", "Hadron", "PID", "Misc"};

  const char *pname[3] = {"PbSc West", "PbSc East", "PbGl"};
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  const double PI = TMath::Pi();
  const double DeltaEta = 0.5;
  //const double NBBC =  3.48e11;  // from DAQ
  const double NBBC =  3.52e11;  // from rejection power
  const double XBBC = 32.51e9;
  const double eXBBC = 3.24e9;

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

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-DC3sigma.root");

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
    double dummy, xpt, xsec[3], exsec[3];  // part
    double rsys[nsys][3], ersys[nsys][3];  // isys-1, part

    for(int part=0; part<3; part++)
    {
      int evtype = 2;
      if(ipt < 22)  // <14GeV use ERT_4x4c
        evtype = 2;
      else  // >14GeV use ERT_4x4b
        evtype = 1;

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

      for(int isys=0; isys<=nsys; isys++)
      {
        double FactorPhoton = 1.;
        double FactorPion = 1.;
        double eFactor = 1e-9;
        double Pile[3] = {1.10, 1.08, 1.05};
        double ePile = 0.02;
        double TrigBBC = 0.91;
        double eTrigBBC = 0.01;
        double Prob = 0.96;
        double eProb = 0.02;
        double ToF[3] = {0.992, 0.992, 0.997};
        double eToF[3] = {0.002, 0.002, 0.002};
        double Conv[3] = {0.849, 0.971, 0.959};
        double eConv[3] = {0.027, 0.023, 0.029};
        double Norm[3] = {0.320, 0.321, 0.250};
        double eNorm[3] = {0.005, 0.007, 0.005};
        double A = 0.28;
        double eA = 0.05;

        double Acc, eAcc, TrigERT, eTrigERT, Miss, eMiss, MissEta, eMissEta,
               Merge1, eMerge1, Merge2, eMerge2, BadPass, eBadPass, Veto, eVeto;
        qt_acc->Query(ipt, part, dummy, Acc, eAcc);
        qt_ert->Query(ipt, part, xpt, TrigERT, eTrigERT);
        qt_miss->Query(ipt, part, dummy, Miss, eMiss);
        qt_miss_eta->Query(ipt, part, dummy, MissEta, eMissEta);
        qt_veto->Query(ipt, part, dummy, Veto, eVeto);

        if(xpt > 10.)
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

        if(xpt > 12.)
        {
          qt_merge1->Query(ipt, part, dummy, Merge1, eMerge1);
          qt_merge2->Query(ipt, part, dummy, Merge2, eMerge2);
          qt_badpass->Query(ipt, part/2, dummy, BadPass, eBadPass);
        }
        else
        {
          Merge1 = eMerge1 = Merge2 = eMerge2 = BadPass = eBadPass = 0.;
        }

        switch(isys)
        {
          case 1:
            double Eff = Conv[part]*Prob*ToF[part];
            double SysPhoton, eSysPhoton;
            qt_sys->Query(ipt<20?ipt/4*4:ipt, 0, dummy, SysPhoton, eSysPhoton);
            FactorPhoton += SysPhoton;
            eFactor = eSysPhoton/Eff;
            break;
          case 2:
            double Eff = Conv[part]*Prob*ToF[part];
            double SysPion, eSysPion;
            qt_sys->Query(ipt<20?ipt/4*4:ipt, 1, dummy, SysPion, eSysPion);
            FactorPion += SysPion;
            eFactor = eSysPion/Eff/Eff;
            break;
          case 3:
            double rbg, erbg;
            qt_rbg->Query(ipt, part/2, dummy, rbg, erbg);
            if(!TMath::Finite(rbg+erbg) || xpt > 10.)
              rbg = erbg = 0.05;
            FactorPion += rbg;
            break;
          case 4:
            A += eA;
            break;
          case 5:
            Prob += eProb;
            break;
          case 6:
            ToF[part] += eToF[part];
          case 7:
            FactorPion += 0.02;
            break;
          case 8:
            Acc += eAcc;
            break;
          case 9:
            Miss += eMiss;
            break;
          case 10:
            MissEta += eMissEta;
            break;
          case 11:
            Veto += eVeto;
            break;
          case 12:
            Merge1 += eMerge1;
            Merge2 += eMerge2;
            break;
          case 13:
            BadPass += eBadPass;
            break;
          case 14:
            Conv[part] += eConv[part];
            break;
          case 15:
            TrigBBC += eTrigBBC;
            break;
          case 16:
            TrigERT += eTrigERT;
            break;
          case 17:
            Norm[part] += eNorm[part];
            break;
          case 18:
            Pile[part] += ePile;
            break;
        }

        double Eff = Conv[part]*Prob*ToF[part];
        double AIso = A*(Veto + MissEta)/(1 + 2*MissEta)*(1 + 2*Miss + Merge1);
        double nbg = (1 + Merge1*Conv[part]*(1 - Conv[part]))*nisoboth + (Miss + AIso)*nisopair + Merge2/2*BadPass*nisopair2pt;
        double e2nbg = pow((1 + Merge1*Conv[part]*(1 - Conv[part]))*enisoboth,2) + pow((Miss + AIso)*enisopair,2) + pow(Merge2/2*BadPass*enisopair2pt,2);

        double ndir = nphoton/Eff*FactorPhoton - nbg/Eff/Eff*FactorPion;
        double endir = sqrt(enphoton*enphoton + e2nbg/Eff/Eff) / Eff;
        if(xpt > 14.)  // >14GeV use ERT_4x4b
        {
          ndir *= Norm[part];
          endir *= Norm[part];
        }
        eFactor = eFactor/ndir;

        qt_pt->Query(ipt, 0, dummy, xpt, dummy);
        double yy = (XBBC/NBBC) / (2*PI*xpt) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
          * ndir / Acc / TrigERT / TrigBBC * Pile[part];
        double eyy = yy*endir/ndir;
        if(isys)
        {
          rsys[isys-1][part] = yy/xsec[part];
          ersys[isys-1][part] = rsys[isys-1][part]*eFactor;
          rsys[isys-1][part] = fabs(rsys[isys-1][part] - 1);
        }
        else if( TMath::Finite(yy+eyy) )
        {
          xsec[part] = yy;
          exsec[part] = eyy;
          qt_cross->Fill(ipt, part, xpt, yy, eyy);
        }
      } // isys
    } // part

    double xbar = xsec[2];
    double exbar = exsec[2];
    if(xpt < 30.)
      Chi2Fit(3, xsec, exsec, xbar, exbar);
    if( TMath::Finite(xbar + exbar) )
      qt_cross->Fill(ipt, 3, xpt, xbar, exbar);
    if( xpt > 6. )
    {
      ndata += xbar;
      nfit += cross_ph->Eval(xpt);
    }

    double rsum = 0.;
    double ersum = 0.;
    for(int igsys=1; igsys<=ngsys; igsys++)
    {
      double rgsys = 0.;
      double ergsys = 0.;
      for(int isys=gsys[igsys-1]; isys<gsys[igsys]; isys++)
      {
        double rbar = TMath::MaxElement(3, rsys[isys-1]);
        double erbar = TMath::MaxElement(3, ersys[isys-1]);
        if(igsys==1)
          Chi2Fit(3, rsys[isys-1], ersys[isys-1], rbar, erbar);
        if( TMath::Finite(rbar+erbar) )
        {
          rgsys += rbar*rbar;
          ergsys += erbar*erbar;
          rsum += rbar*rbar;
          ersum += erbar*erbar;
        }
      } // isys
      if( TMath::Finite(rgsys+ergsys) )
        qt_cross->Fill(ipt, 4+igsys, xpt, sqrt(rgsys), sqrt(ergsys));
    } // igsys
    if( TMath::Finite(rsum+ersum) )
      qt_cross->Fill(ipt, 4, xpt, sqrt(rsum), sqrt(ersum));
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
    aset(gr, "p_{T} [GeV]", "Ed^{3}#sigma/dp^{3} [pb GeV^{-2} c^{3}]", 6.1,30., 1e-1, 2e3);
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
  legi(1, 0.2,0.8,0.7,0.9);
  leg1->SetNColumns(3);
  for(int igsys=0; igsys<=ngsys; igsys++)
  {
    TGraphErrors *gr = qt_cross->Graph(4+igsys);
    aset(gr, "p_{T} [GeV]", "SysErr", 6.1,30., 0.,0.2);
    style(gr, igsys+20, igsys+1);
    gr->SetLineStyle(igsys/3*8+1);
    char *opt = igsys==0 ? "AP" : "L";
    gr->Draw(opt);
    leg1->AddEntry(gr, sysname[igsys], igsys==0?"P":"L");
  }
  leg1->Draw();
  c10->Print("plots/SysErr-isophoton.pdf");

  qt_cross->Write();
  for(int part=0; part<3; part++)
  {
    mcw( part, Form("Minv-isoboth-part%d",part) );
    mcw( part+3, Form("Minv-isopair-part%d",part) );
    mcw( part+6, Form("Minv-isopair2pt-part%d",part) );
  }
  qt_cross->Close();
}
