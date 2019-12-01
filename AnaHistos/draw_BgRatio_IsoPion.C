#include "GlobalVars.h"
#include "QueryTree.h"
#include "BgGPR.h"

void draw_Each(const QueryTree *qt_rbg, const int beam, const int isotype, const int pttype)
{
  const double Prob = 0.96;
  const double eProb = 0.02;
  const double ToF = 0.995;
  const double eToF = 0.005;
  const double Conv = 0.90;
  const double eConv = 0.05;
  const double A = 0.28;
  const double eA = 0.05;

  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");
  QueryTree *qt_miss_eta = new QueryTree("data/MissingRatio-eta.root");
  QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root");
  QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root");
  QueryTree *qt_badpass = new QueryTree("data/MergePassRate.root");
  QueryTree *qt_veto = new QueryTree("data/SelfVeto.root");

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
    double xpt, yy[3], eyy[3];

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

      double xpt, Miss, eMiss, MissEta, eMissEta, Merge1, eMerge1, Merge2, eMerge2, BadPass, eBadPass, Veto, eVeto;
      qt_miss->Query(ipt, 3, xpt, Miss, eMiss);
      qt_miss_eta->Query(ipt, 3, xpt, MissEta, eMissEta);
      qt_merge1->Query(ipt, 3, xpt, Merge1, eMerge1);
      qt_merge2->Query(ipt, 3, xpt, Merge2, eMerge2);
      qt_badpass->Query(ipt, 2, xpt, BadPass, eBadPass);
      qt_veto->Query(ipt, 3, xpt, Veto, eVeto);

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

      double Eff = Conv[part] * Prob * ToF[part];
      double AIso = A * (Veto+MissEta)/(1.+2.*MissEta) * (1+2.*Miss+Merge1);
      double ndir = nphoton/Eff - (1. + Merge1*Conv[part]*(1.-Conv[part])) * nisoboth/Eff/Eff - Miss * nisopair/Eff/Eff - Merge2/2.*BadPass * nisopair2pt/Eff/Eff - AIso * nisopair/Eff/Eff;
      double endir = sqrt(pow(enphoton,2)/pow(Conv[part],2) + (pow(enisoboth,2)*pow(1. + (1. - Conv[part])*Conv[part]*Merge1,2))/pow(Conv[part],4) + 0.25*pow(BadPass,2)*pow(enisopair2pt,2)*pow(Merge2,2) + (pow(A,2)*pow(eVeto,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(nisopair,2))/pow(1. + 2.*Miss,2) + 0.25*pow(BadPass,2)*pow(eMerge2,2)*pow(nisopair2pt,2) + 0.25*pow(eBadPass,2)*pow(Merge2,2)*pow(nisopair2pt,2) + pow(eConv[part],2)*pow(-((((1. - Conv[part])*Merge1 - Conv[part]*Merge1)*nisoboth)/pow(Conv[part],2)) + (2*(1. + (1. - Conv[part])*Conv[part]*Merge1)*nisoboth)/pow(Conv[part],3) + (2*Miss*nisopair)/pow(Conv[part],3) - nphoton/pow(Conv[part],2),2) + (pow(eA,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(nisopair,2)*pow(Miss + Veto,2))/pow(1. + 2.*Miss,2) + pow(enisopair,2)*pow(-(Miss/pow(Conv[part],2)) - (A*(1 + Merge1 + 2.*Miss)*(Miss + Veto))/(1. + 2.*Miss),2) + pow(eMerge1,2)*pow(-(((1. - Conv[part])*nisoboth)/Conv[part]) - (A*nisopair*(Miss + Veto))/(1. + 2.*Miss),2) + pow(eMiss,2)*pow(-(nisopair/pow(Conv[part],2)) - (A*(1 + Merge1 + 2.*Miss)*nisopair)/(1. + 2.*Miss) - (2.*A*nisopair*(Miss + Veto))/(1. + 2.*Miss) + (2.*A*(1 + Merge1 + 2.*Miss)*nisopair*(Miss + Veto))/pow(1. + 2.*Miss,2),2));

      if(ipt >= 22)  // >14GeV use ERT_4x4b
      {
        ndir *= Norm[part];
        endir *= Norm[part];
      }

      yy[part] = (XBBC/NBBC) / (2*PI*xpt) / (pTbin[ipt+1]-pTbin[ipt]) / DeltaEta
        * ndir / Acc / TrigERT / TrigBBC * Pile[part];
      eyy[part] = yy[part] * sqrt( pow(endir/ndir,2)
          + pow(eAcc/Acc,2)
          + pow(eTrigERT/TrigERT,2)
          + pow(eToF[part]/ToF[part],2)
          + pow(eProb/Prob,2)
          + pow(ePile/Pile[part],2)
          //+ pow(eTrigBBC/TrigBBC,2) + pow(eXBBC/XBBC,2)
          );
      if( TMath::Finite(yy[part]+eyy[part]) )
        qt_cross->Fill(ipt, part, xpt, yy[part], eyy[part]);
    } // part

    double ybar, eybar;
    Chi2Fit(3, yy, eyy, ybar, eybar);
    if( TMath::Finite(ybar + eybar) )
      qt_cross->Fill(ipt, 3, xpt, ybar, eybar);
    if( ipt >= 12 )
    {
      ndata += ybar;
      nfit += cross_ph->Eval(xpt);
    }
  } // ipt
}

void draw_Each_0(const QueryTree *qt_rbg, const int beam, const int isotype, const int pttype)
{
  const int checkmap = 1;
  const char *ptname = pttype ? "2pt" : "";

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  TH2 *h2_pion = (TH2*)f->Get("h2_2photon_pol_0");
  h2_pion = (TH2*)h2_pion->Clone();
  h2_pion->Reset();

  const unsigned nData = 45;
  vector<double> x(nData), y(nData), sigma_y(nData);

  for(int icr=0; icr<2; icr++)
  {
    h2_pion->Reset();
    for(int iso=0; iso<2; iso++)
      for(int ipol=0; ipol<2; ipol++)
      {
        int ih = beam + 3*icr + 3*2*ipol + 3*2*2*checkmap + 3*2*2*2*(1-isotype*iso) + 3*2*2*2*2*(1-(1-isotype)*iso);
        TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon%s_pol_%d",ptname,ih));
        h2_pion->Add(h2_tmp);
      }

    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      int ptbin_first = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt]);
      int ptbin_last = h2_pion->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;
      TH1 *h_minv = h2_pion->ProjectionY("h_minv", ptbin_first,ptbin_last);

      x.clear();
      y.clear();
      sigma_y.clear();

      for(int biny=68; biny<=87; biny++)
      {
        double xx = h_minv->GetXaxis()->GetBinCenter(biny);
        double yy = h_minv->GetBinContent(biny);
        double sigma_yy = h_minv->GetBinError(biny);
        x.push_back(xx);
        y.push_back(yy);
        sigma_y.push_back(sigma_yy);
      }

      for(int biny=188; biny<=212; biny++)
      {
        double xx = h_minv->GetXaxis()->GetBinCenter(biny);
        double yy = h_minv->GetBinContent(biny);
        double sigma_yy = h_minv->GetBinError(biny);
        x.push_back(xx);
        y.push_back(yy);
        sigma_y.push_back(sigma_yy);
      }

      double nbg, dnbg;
      BgGPR(x, y, sigma_y, nbg, dnbg);
      nbg /= 0.001;
      dnbg /= 0.001;

      double npeak = h_minv->Integral(113, 162);
      double rbg = nbg/npeak;
      double erbg = rbg*sqrt(dnbg*dnbg/nbg/nbg + 1./npeak);
      if( !TMath::Finite(rbg+erbg) || rbg < 0. || erbg > 1. )
      {
        rbg = 0.;
        erbg = 1.;
      }

      double xpt = (pTbin_pol[ipt] + pTbin_pol[ipt+1]) / 2.;
      int part = beam + 3*isotype + 3*2*pttype + 3*2*2*icr;
      qt_rbg->Fill(ipt, part, xpt, rbg, erbg);
    } // ipt
  } // icr
}

void draw_BgRatio_IsoPion()
{
  gSystem->Load("libGausProc.so");

  QueryTree *qt_rbg = new QueryTree("data/BgRatio-isophoton.root", "RECREATE");

  for(int beam=0; beam<3; beam++)
    for(int isotype=0; isotype<2; isotype++)
      for(int pttype=0; pttype<2; pttype++)
        draw_Each(qt_rbg, beam, isotype, pttype);

  qt_rbg->Save();
}
