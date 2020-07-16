#include "GlobalVars.h"
#include "QueryTree.h"
#include "Ftest.h"
#include "Chi2Fit.h"
#include "IsoPhotonALL.h"

void draw_IsoPhotonALL()
{
  const char *beam_list[3] = {"A_{L}^{Blue}", "A_{L}^{Yellow}", "A_{LL}"};
  const char *region[3] = {"Sig+BG", "Sig", "Combined"};
  const char *crossing_list[2] = {"Even", "Odd"};
  const char *pattern_list[4] = {"SOOSSOO", "OSSOOSS", "SSOO", "OOSS"};

  QueryTree *qt_all = new QueryTree("data/IsoPhotonALL.root", "RECREATE");

  QueryTree *qt_asym = new QueryTree("data/isophoton-asym-fill-tightcut.root");
  qt_asym->SetQuiet();
  int imul = 1;

  QueryTree *qt_allpion = new QueryTree("data/IsoPionALL.root");
  QueryTree *qt_rbg = new QueryTree("data/BgRatio-isophoton.root");

  vector<double> *vp_ALL = new vector<double>[8];
  vector<double> *vp_eALL = new vector<double>[8];

  cout.precision(4);
  for(int beam=0; beam<3; beam++)
  {
    cout << "beam " << beam << endl;

    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      for(int icr=0; icr<2; icr++)
        for(int pattern=0; pattern<4; pattern++)
        {
          int id = icr + 2*pattern;
          vp_ALL[id].clear();
          vp_eALL[id].clear();

          int ig = imul + 6*beam + 6*3*icr + 6*3*2*pattern + 6*3*2*4*ipt;
          TSQLResult *res = qt_asym->Query(ig); // runnumber:runnumber:value:error:errorlow:errorhigh
          TSQLRow *row;
          while( row = res->Next() )
          {
            TString field2 = row->GetField(2);
            TString field3 = row->GetField(3);
            double ALL = field2.Atof();
            double eALL = field3.Atof();

            vp_ALL[id].push_back(ALL);
            vp_eALL[id].push_back(eALL);
            delete row;
          }
          delete res;
        } // icr, pattern

      double F, p;
      Ftest(8, vp_ALL, vp_eALL, F, p);
      cout << pTbin_pol[ipt] << "-" << pTbin_pol[ipt+1] << " & " << F << " & " << p << " \\\\" << endl;
    } // ipt
  } // beam

  TF1 *fn_mean = new TF1("fn_mean", "pol0");

  for(int beam=0; beam<3; beam++)
    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      double sig[2][8] = {};
      double esig[2][8] = {};
      for(int icr=0; icr<2; icr++)
        for(int pattern=0; pattern<4; pattern++)
        {
          double xpt, allpion[2], eallpion[2], rbg[3], erbg[3];
          for(int pttype=0; pttype<2; pttype++)
          {
            int part = beam + 3*pttype + ngr_pion/2*3;
            qt_allpion->Query(ipt, part, xpt, allpion[pttype], eallpion[pttype]);
            if( !TMath::Finite(allpion[pttype]) )
              allpion[pttype] = 0.;
            if( !TMath::Finite(eallpion[pttype]) )
              eallpion[pttype] = 0.;
          }
          for(int ibg=0; ibg<3; ibg++)
          {
            part = ibg + 3*icr;
            qt_rbg->Query(ipt, part, xpt, rbg[ibg], erbg[ibg]);
            if( !TMath::Finite(rbg[ibg]) )
              rbg[ibg] = 0.;
            if( !TMath::Finite(erbg[ibg]) )
              erbg[ibg] = 0.;
          }

          int igr = beam + 3*icr + 3*2*pattern;
          if(ipt == 0)
            mc(igr, 4,4);
          mcd(igr, ipt+1);

          int ig = imul + 6*beam + 6*3*icr + 6*3*2*pattern + 6*3*2*4*ipt;
          TGraphErrors *gr = qt_asym->Graph(ig); 
          gr->SetTitle(Form("p_{T}: %.1f-%.1f GeV",pTbin_pol[ipt],pTbin_pol[ipt+1]));
          aset(gr, "runnumber",beam_list[beam], 386700.,398200., -0.5,0.5);
          style(gr, 20, 1);
          gr->Draw("AP");

          gr->Fit(fn_mean, "Q");
          double mean = fn_mean->GetParameter(0);
          double emean = fn_mean->GetParError(0);
          if( TMath::Finite(mean+emean) )
            qt_all->Fill(ipt, igr, xpt, mean, emean);

          for(int isys=0; isys<2; isys++)
          {
            int ipat = icr + 2*pattern;
            for(int ibg=0; ibg<3; ibg++)
              rbg[ibg] *= (1 + 0.04*isys);
            sig[isys][ipat] = (mean - rbg[0]*allpion[0] - rbg[1]*allpion[1]) / (1 - rbg[0] - rbg[1] - rbg[2]);
            esig[isys][ipat] = sqrt( emean*emean + pow(rbg[0]*eallpion[0],2) + pow(rbg[1]*eallpion[1],2) ) / (1 - rbg[0] - rbg[1] - rbg[2]);

            int igr = beam + 3*icr + 3*2*pattern + ngr_photon + ngr_photon*3*isys;
            qt_all->Fill(ipt, igr, xpt, sig[isys][ipat], esig[isys][ipat]);
          } // isys
        } // icr, pattern

      for(int isys=0; isys<2; isys++)
      {
        double comb, ecomb;
        Chi2Fit(8, sig[isys], esig[isys], comb, ecomb);
        if( TMath::Finite(comb+ecomb) )
        {
          int igr = beam + ngr_photon*2 + ngr_photon*3*isys;
          qt_all->Fill(ipt, igr, xpt, comb, ecomb);
        }
      } // isys
    } // beam, ipt

  for(int beam=0; beam<3; beam++)
  {
    int ic = beam + ngr_photon;
    mc(ic, 2,2);
    legi(0, 0.2,0.8,0.7,0.9);
    leg0->SetNColumns(2);
    leg0->SetTextSize(0.02);

    for(int ibg=0; ibg<2; ibg++)
      for(int icr=0; icr<2; icr++)
        for(int pattern=0; pattern<4; pattern++)
        {
          mcd(ic, ibg+1);
          int igr = beam + 3*icr + 3*2*pattern + 3*2*4*ibg;
          TGraphErrors *gr_all = qt_all->Graph(igr);

          gr_all->SetTitle( Form("#gamma^{dir} %s %s",beam_list[beam],region[ibg]) );
          aset(gr_all, "p_{T} [GeV]",beam_list[beam], 0.,20., -0.2,0.4);
          style(gr_all, 1+icr, 1+pattern);
          gr_all->SetMarkerSize(0.);
          if(icr==0 && pattern==0)
            gr_all->Draw("AP");
          else
            gr_all->Draw("P");
          if(ibg == 0)
            leg0->AddEntry(gr_all, Form("%s %s",pattern_list[pattern],crossing_list[icr]), "L");
        } // ibg, icr, pattern
    leg0->Draw();

    mcd(ic, 3);
    gPad->SetGridy();
    int igr = beam + ngr_photon*2;
    TGraphErrors *gr_all = qt_all->Graph(igr);
    gr_all->SetTitle( Form("#gamma^{dir} %s %s",beam_list[beam],region[2]) );
    aset(gr_all, "p_{T} [GeV]",beam_list[beam], 0.,20., -0.06,0.05);
    style(gr_all, 1, 1);
    gr_all->SetMarkerSize(0.);
    gr_all->Draw("AP");

    qt_all->Write();
    for(int icr=0; icr<2; icr++)
      for(int pattern=0; pattern<4; pattern++)
      {
        int igr = beam + 3*icr + 3*2*pattern;
        mcw( igr, Form("beam%d-pattern%d-cross%d",
              beam, pattern, icr) );
      } // icr, pattern
    int igr = beam + ngr_photon;
    mcw( igr, Form("beam%d-combined",
          beam) );
  } // beam

  qt_all->Close();
}
