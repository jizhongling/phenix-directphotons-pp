#include "GlobalVars.h"
#include "QueryTree.h"
#include "Ftest.h"
#include "Chi2Fit.h"
#include "IsoPhotonALL.h"

void draw_Pi0ALL()
{
  const char *beam_list[3] = {"A_{L}^{Blue}", "A_{L}^{Yellow}", "A_{LL}"};
  const char *region[4] = {"Sig+BG", "BG", "Sig", "Combined"};
  const char *crossing_list[2] = {"Even", "Odd"};
  const char *pattern_list[4] = {"SOOSSOO", "OSSOOSS", "SSOO", "OOSS"};

  QueryTree *qt_all = new QueryTree("data/Pi0ALL.root", "RECREATE");

  TGraphErrors *gr_inseok[3];
  gr_inseok[0] = new TGraphErrors("data/inseok-alblue.txt", "%lg %lg %lg");
  gr_inseok[1] = new TGraphErrors("data/inseok-alyellow.txt", "%lg %lg %lg");
  gr_inseok[2] = new TGraphErrors("data/inseok-all.txt", "%lg %lg %lg");

  QueryTree *qt_asym = new QueryTree("data/pion-asym.root");
  qt_asym->SetQuiet();

  QueryTree *qt_rbg = new QueryTree("data/BgRatio-pion.root");

  vector<double> *vp_ALL = new vector<double>[8];
  vector<double> *vp_eALL = new vector<double>[8];

  cout.precision(4);
  for(int beam=0; beam<3; beam++)
    for(int ibg=0; ibg<2; ibg++)
    {
      cout << "beam " << beam << ", bg " << ibg << endl;
      for(int ipt=0; ipt<npT_pol; ipt++)
      {
        for(int icr=0; icr<2; icr++)
          for(int pattern=0; pattern<4; pattern++)
          {
            int id = icr + 2*pattern;
            vp_ALL[id].clear();
            vp_eALL[id].clear();

            int ig = beam + 3*ibg + 3*2*icr + 3*2*2*pattern + 3*2*2*4*ipt;
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
    } // beam, ibg

  TF1 *fn_mean = new TF1("fn_mean", "pol0");

  for(int beam=0; beam<3; beam++)
    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      double sig[8] = {};
      double esig[8] = {};
      for(int icr=0; icr<2; icr++)
        for(int pattern=0; pattern<4; pattern++)
        {
          double xpt, rbg, erbg;
          qt_rbg->Query(ipt, icr, xpt, rbg, erbg);
          double mean[2], emean[2];

          for(int ibg=0; ibg<2; ibg++)
          {
            int igr = beam + 3*icr + 3*2*pattern + 3*2*4*ibg;
            if(ipt == 0)
              mc(igr, 4,4);
            mcd(igr, ipt+1);

            int ig = beam + 3*ibg + 3*2*icr + 3*2*2*pattern + 3*2*2*4*ipt;
            TGraphErrors *gr = qt_asym->Graph(ig); 
            gr->SetTitle(Form("p_{T}: %.1f-%.1f GeV",pTbin_pol[ipt],pTbin_pol[ipt+1]));
            aset(gr, "runnumber","A_{LL}", 386700.,398200., -0.5,0.5);
            style(gr, 20, 1);
            gr->Draw("AP");

            gr->Fit(fn_mean, "Q");
            mean[ibg] = fn_mean->GetParameter(0);
            emean[ibg] = fn_mean->GetParError(0);
            if( TMath::Finite(mean[ibg]+emean[ibg]) )
              qt_all->Fill(ipt, igr, xpt, mean[ibg], emean[ibg]);
          } // ibg

          int ipat = icr + 2*pattern;
          sig[ipat] = (mean[0] - rbg*mean[1])/(1 - rbg);
          esig[ipat] = sqrt(pow(emean[0],2)/pow(1 - rbg,2) + (pow(emean[1],2)*pow(rbg,2))/pow(1 - rbg,2) + pow(erbg,2)*pow(-(mean[1]/(1 - rbg)) + (mean[0] - mean[1]*rbg)/pow(1 - rbg,2),2));

          int igr = beam + 3*icr + 3*2*pattern + ngr_pi0;
          qt_all->Fill(ipt, igr, xpt, sig[ipat], esig[ipat]);
        } // icr, pattern

      double comb, ecomb;
      Chi2Fit(8, sig, esig, comb, ecomb);
      if( TMath::Finite(comb+ecomb) )
      {
        int igr = beam + ngr_pi0/2*3;
        qt_all->Fill(ipt, igr, xpt, comb, ecomb);
      }
    } // beam, ipt

  for(int beam=0; beam<3; beam++)
  {
    int ic = beam + ngr_pi0;
    mc(ic, 2,2);
    legi(0, 0.2,0.8,0.7,0.9);
    leg0->SetNColumns(2);
    leg0->SetTextSize(0.02);

    for(int ibg=0; ibg<3; ibg++)
      for(int icr=0; icr<2; icr++)
        for(int pattern=0; pattern<4; pattern++)
        {
          mcd(ic, ibg+1);
          int igr = beam + 3*icr + 3*2*pattern + 3*2*4*ibg;
          TGraphErrors *gr_all = qt_all->Graph(igr);

          gr_all->SetTitle( Form("#pi^{0} %s %s",beam_list[beam],region[ibg]) );
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

    mcd(ic, 4);
    legi(1, 0.2,0.8,0.6,0.9);
    gPad->SetGridy();
    int igr = beam + ngr_pi0/2*3;
    TGraphErrors *gr_all = qt_all->Graph(igr);
    gr_all->SetTitle( Form("#pi^{0} %s %s",beam_list[beam],region[3]) );
    aset(gr_all, "p_{T} [GeV]",beam_list[beam], 0.,20., -0.02,0.03);
    style(gr_all, 1, 1);
    style(gr_inseok[beam], 1, 2);
    gr_all->SetMarkerSize(0.);
    gr_inseok[beam]->SetMarkerSize(0.);
    gr_all->Draw("AP");
    gr_inseok[beam]->Draw("P");
    leg1->AddEntry(gr_all, "My check", "L");
    leg1->AddEntry(gr_inseok[beam], "Inseok", "L");
    leg1->Draw();

    qt_all->Write();
    for(int ibg=0; ibg<2; ibg++)
      for(int icr=0; icr<2; icr++)
        for(int pattern=0; pattern<4; pattern++)
        {
          int igr = beam + 3*icr + 3*2*pattern + 3*2*4*ibg;
          mcw( igr, Form("beam%d-bg%d-pattern%d-cross%d", beam, ibg, pattern, icr) );
        } // ibg, icr, pattern
    int igr = beam + ngr_pi0;
    mcw( igr, Form("beam%d-combined", beam) );
  } // beam

  qt_all->Close();
}
