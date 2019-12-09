#include "GlobalVars.h"
#include "QueryTree.h"
#include "FitMinv.h"

void draw_BgRatio_IsoPhoton()
{
  const double meff[2] = {0.98, 0.99};
  const double bck = 1.15;
  const double Prob = 0.96;
  const double eProb = 0.02;
  const double ToF = 0.995;
  const double eToF = 0.005;
  const double Conv = 0.90;
  const double eConv = 0.05;
  const double A = 0.28;
  const double eA = 0.05;

  QueryTree *qt_rbg = new QueryTree("data/BgRatio-isophoton.root", "RECREATE");

  QueryTree *qt_miss = new QueryTree("data/MissingRatio.root");
  QueryTree *qt_miss_eta = new QueryTree("data/MissingRatio-eta.root");
  QueryTree *qt_merge1 = new QueryTree("data/Merge-1photon.root");
  QueryTree *qt_merge2 = new QueryTree("data/Merge-2photon.root");
  QueryTree *qt_badpass = new QueryTree("data/MergePassRate.root");
  QueryTree *qt_veto = new QueryTree("data/SelfVeto.root");

  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  TH1 *h_photon = (TH1*)f->Get("h_1photon_pol_0");
  h_photon = (TH1*)h_photon->Clone();

  TH2 *h2_pion[2];
  h2_pion[0] = (TH2*)f->Get("h2_2photon_pol_0");
  h2_pion[0] = (TH2*)h2_pion[0]->Clone();
  h2_pion[1] = (TH2*)h2_pion[0]->Clone();

  int beam = 2;
  int checkmap = 1;

  for(int icr=0; icr<2; icr++)
  {
    h_photon->Reset();
    for(int ipol=0; ipol<2; ipol++)
    {
      int ih = beam + 3*icr + 3*2*ipol + 3*2*2*checkmap;
      TH1 *h_tmp = (TH1*)f->Get(Form("h_1photon_pol_%d",ih));
      h_photon->Add(h_tmp);
    } // ipol

    for(int pttype=0; pttype<2; pttype++)
    {
      const char *ptname = pttype ? "2pt" : "";
      for(int isotype=0; isotype<2; isotype++)
      {
        int ic = isotype + 2*pttype + 2*2*icr;
        mc(ic, 4,4);
        h2_pion[isotype]->Reset();
        for(int iso=0; iso<2; iso++)
          for(int ipol=0; ipol<2; ipol++)
          {
            int ih = beam + 3*icr + 3*2*ipol + 3*2*2*checkmap + 3*2*2*2*(1-isotype*iso) + 3*2*2*2*2*(1-(1-isotype)*iso);
            TH2 *h2_tmp = (TH2*)f->Get(Form("h2_2photon%s_pol_%d",ptname,ih));
            h2_pion[isotype]->Add(h2_tmp);
          } // ipol
      } // isotype

      for(int ipt=0; ipt<npT_pol; ipt++)
      {
        double xpt = (pTbin_pol[ipt] + pTbin_pol[ipt+1]) / 2.;
        int ptbin_first = h_photon->GetXaxis()->FindBin(pTbin_pol[ipt]);
        int ptbin_last = h_photon->GetXaxis()->FindBin(pTbin_pol[ipt+1]) - 1;

        double nphoton = h_photon->Integral(ptbin_first,ptbin_last);
        double enphoton = sqrt(nphoton);

        double npion[2], enpion[2];
        for(int isotype=0; isotype<2; isotype++)
        {
          TH1 *h_minv = h2_pion[isotype]->ProjectionY("h_minv", ptbin_first,ptbin_last);

          int ic = isotype + 2*pttype + 2*2*icr;
          mcd(ic, ipt+1);
          h_minv->Rebin(10);
          h_minv->SetTitle( Form("p_{T}: %3.1f-%3.1f GeV", pTbin_pol[ipt], pTbin_pol[ipt+1]) );
          if(xpt < 10.)  // <10GeV +-25MeV; >10GeV +-35MeV
          {
            FitMinv(h_minv, npion[isotype], enpion[isotype], true, 0.11,0.16);
            npion[isotype] /= meff[0];
            enpion[isotype] /= meff[0];
          }
          else if(xpt < 16.)  // <16GeV subtract background
          {
            FitMinv(h_minv, npion[isotype], enpion[isotype], true, 0.10,0.17);
            npion[isotype] /= meff[1];
            enpion[isotype] /= meff[1];
          }
          else  // >16GeV don't subtract background
          {
            FitMinv(h_minv, npion[isotype], enpion[isotype], false, 0.10,0.17);
            npion[isotype] /= meff[1] * bck;
            enpion[isotype] /= meff[1] * bck;
          }
          delete h_minv;
        } // isotype

        double Miss, eMiss, MissEta, eMissEta, Merge1, eMerge1, Merge2, eMerge2, BadPass, eBadPass, Veto, eVeto;
        qt_miss->Query(ipt, 3, xpt, Miss, eMiss);
        qt_miss_eta->Query(ipt, 3, xpt, MissEta, eMissEta);
        qt_merge1->Query(ipt, 3, xpt, Merge1, eMerge1);
        qt_merge2->Query(ipt, 3, xpt, Merge2, eMerge2);
        qt_badpass->Query(ipt, 2, xpt, BadPass, eBadPass);
        qt_veto->Query(ipt, 3, xpt, Veto, eVeto);
        if(xpt < 12.)
        {
          Merge1 = 0.;
          eMerge1 = 0.;
          Merge2 = 0.;
          eMerge2 = 0.;
          BadPass = 0.;
          eBadPass = 0.;
        }

        int part = 3*icr;
        if(pttype == 0)
        {
          double Eff = Conv*Prob*ToF;

          double rbg = (1 + Merge1*(1 - Conv)/Eff)*npion[0]/nphoton + Miss/Eff*npion[1]/nphoton;
          double erbg = sqrt(pow(eToF,2)*pow(-(((1 - Conv)*Merge1*npion[0])/(Conv*nphoton*Prob*pow(ToF,2))) - (Miss*npion[1])/(Conv*nphoton*Prob*pow(ToF,2)),2) + pow(eProb,2)*pow(-(((1 - Conv)*Merge1*npion[0])/(Conv*nphoton*pow(Prob,2)*ToF)) - (Miss*npion[1])/(Conv*nphoton*pow(Prob,2)*ToF),2) + (pow(enpion[0],2)*pow(1 + ((1 - Conv)*Merge1)/(Conv*Prob*ToF),2))/pow(nphoton,2) + pow(enphoton,2)*pow(-((npion[0]*(1 + ((1 - Conv)*Merge1)/(Conv*Prob*ToF)))/pow(nphoton,2)) - (Miss*npion[1])/(Conv*pow(nphoton,2)*Prob*ToF),2) + pow(eConv,2)*pow((npion[0]*(-(((1 - Conv)*Merge1)/(pow(Conv,2)*Prob*ToF)) - Merge1/(Conv*Prob*ToF)))/nphoton - (Miss*npion[1])/(pow(Conv,2)*nphoton*Prob*ToF),2) + (pow(enpion[1],2)*pow(Miss,2))/(pow(Conv,2)*pow(nphoton,2)*pow(Prob,2)*pow(ToF,2)) + (pow(1 - Conv,2)*pow(eMerge1,2)*pow(npion[0],2))/(pow(Conv,2)*pow(nphoton,2)*pow(Prob,2)*pow(ToF,2)) + (pow(eMiss,2)*pow(npion[1],2))/(pow(Conv,2)*pow(nphoton,2)*pow(Prob,2)*pow(ToF,2)));
          qt_rbg->Fill(ipt, part, xpt, rbg, erbg);

          rbg = A*(1 + 2.*Miss + Merge1)/(1 + 2.*MissEta)*(Veto + MissEta/Eff)*npion[0]/nphoton;
          erbg = sqrt((pow(A,2)*pow(eVeto,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(npion[0],2))/(pow(1 + 2.*MissEta,2)*pow(nphoton,2))+ (pow(A,2)*pow(eToF,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(MissEta,2)*pow(npion[0],2))/(pow(Conv,2)*pow(1 + 2.*MissEta,2)*pow(nphoton,2)*pow(Prob,2)*pow(ToF,4)) + (pow(A,2)*pow(eProb,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(MissEta,2)*pow(npion[0],2))/(pow(Conv,2)*pow(1 + 2.*MissEta,2)*pow(nphoton,2)*pow(Prob,4)*pow(ToF,2)) + (pow(A,2)*pow(eConv,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(MissEta,2)*pow(npion[0],2))/(pow(Conv,4)*pow(1 + 2.*MissEta,2)*pow(nphoton,2)*pow(Prob,2)*pow(ToF,2)) + (pow(A,2)*pow(enpion[0],2)*pow(1 + Merge1 + 2.*Miss,2)*pow(MissEta/(Conv*Prob*ToF) + Veto,2))/(pow(1 + 2.*MissEta,2)*pow(nphoton,2))+ (pow(A,2)*pow(enphoton,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(npion[0],2)*pow(MissEta/(Conv*Prob*ToF) + Veto,2))/(pow(1 + 2.*MissEta,2)*pow(nphoton,4))+ (pow(A,2)*pow(eMerge1,2)*pow(npion[0],2)*pow(MissEta/(Conv*Prob*ToF) + Veto,2))/(pow(1 + 2.*MissEta,2)*pow(nphoton,2))+ (4.*pow(A,2)*pow(eMiss,2)*pow(npion[0],2)*pow(MissEta/(Conv*Prob*ToF) + Veto,2))/(pow(1 + 2.*MissEta,2)*pow(nphoton,2))+ (pow(eA,2)*pow(1 + Merge1 + 2.*Miss,2)*pow(npion[0],2)*pow(MissEta/(Conv*Prob*ToF) + Veto,2))/(pow(1 + 2.*MissEta,2)*pow(nphoton,2))+ pow(eMissEta,2)*pow((A*(1 + Merge1 + 2.*Miss)*npion[0])/(Conv*(1 + 2.*MissEta)*nphoton*Prob*ToF)- (2.*A*(1 + Merge1 + 2.*Miss)*npion[0]*(MissEta/(Conv*Prob*ToF) + Veto))/(pow(1 + 2.*MissEta,2)*nphoton),2));
          qt_rbg->Fill(ipt, part+2, xpt, rbg, erbg);
        }
        else
        {
          double rbg = Merge2*BadPass/2./Prob/Prob*npion[1]/nphoton;
          double erbg = sqrt((1.*pow(BadPass,2)*pow(eProb,2)*pow(Merge2,2)*pow(npion[1],2))/(pow(nphoton,2)*pow(Prob,6)) + (0.25*pow(BadPass,2)*pow(enpion[1],2)*pow(Merge2,2))/(pow(nphoton,2)*pow(Prob,4)) + (0.25*pow(BadPass,2)*pow(enphoton,2)*pow(Merge2,2)*pow(npion[1],2))/(pow(nphoton,4)*pow(Prob,4)) + (0.25*pow(BadPass,2)*pow(eMerge2,2)*pow(npion[1],2))/(pow(nphoton,2)*pow(Prob,4)) + (0.25*pow(eBadPass,2)*pow(Merge2,2)*pow(npion[1],2))/(pow(nphoton,2)*pow(Prob,4)));
          qt_rbg->Fill(ipt, part+1, xpt, rbg, erbg);
        }
      } // ipt
    } // pttype
  } // icr

  qt_rbg->Save();
}
