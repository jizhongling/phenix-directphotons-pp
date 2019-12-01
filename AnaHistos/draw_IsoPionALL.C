#include "GlobalVars.h"
#include "QueryTree.h"
#include "Ftest.h"
#include "Chi2Fit.h"

void draw_Each(const QueryTree *qt_all, const int beam, const int isotype, const int pttype)
{
  static int ncalled = 0;
  const int nstart = 3*2*2*ncalled;
  ncalled++;

  const char *region[4] = {"Sig+BG", "BG", "Sig", "Combined"};
  const char *crossing_list[2] = {"Even", "Odd"};
  const char *pattern_list[4] = {"SOOSSOO", "OSSOOSS", "SSOO", "OOSS"};

  QueryTree *qt_asym = new QueryTree("data/isophoton-asym.root");
  qt_asym->SetQuiet();

  QueryTree *qt_rbg = new QueryTree("data/BgRatio-isophoton.root");

  vector<double> *vp_ALL = new vector<double>[8];
  vector<double> *vp_eALL = new vector<double>[8];

  for(int ibg=0; ibg<2; ibg++)
  {
    cout << "beam " << beam << ", isotype " << isotype << ", pttype " << pttype << ", bg " << ibg << endl;

    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      for(int icr=0; icr<2; icr++)
        for(int pattern=0; pattern<4; pattern++)
        {
          int id = icr + 2*pattern;
          vp_ALL[id].clear();
          vp_eALL[id].clear();

          int imul = 1 + isotype + 2*ibg + 2*2*pttype;
          int ig = imul + 9*beam + 9*3*icr + 9*3*2*pattern + 9*3*2*4*ipt;
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
  } // ibg

  TF1 *fn_mean = new TF1("fn_mean", "pol0");

  for(int ipt=0; ipt<npT_pol; ipt++)
  {
    double sig[8] = {};
    double esig[8] = {};
    for(int icr=0; icr<2; icr++)
      for(int pattern=0; pattern<4; pattern++)
      {
        double xpt, rbg, erbg;
        int part = beam + 3*isotype + 3*2*pttype + 3*2*2*icr;
        qt_rbg->Query(ipt, part, xpt, rbg, erbg);
        double mean[2], emean[2];

        for(int ibg=0; ibg<2; ibg++)
        {
          int igr = icr + 2*pattern + 2*4*ibg + nstart;
          if(ipt == 0)
            mc(igr, 4,4);
          mcd(igr, ipt+1);

          int imul = 1 + isotype + 2*ibg + 2*2*pttype;
          int ig = imul + 9*beam + 9*3*icr + 9*3*2*pattern + 9*3*2*4*ipt;
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
        esig[ipat] = sqrt(emean[0]*emean[0] + rbg*rbg*emean[1]*emean[1])/(1 - rbg);

        int igr = icr + 2*pattern + 2*4*2 + nstart;
        qt_all->Fill(ipt, igr, xpt, sig[ipat], esig[ipat]);
      } // icr, pattern

    double comb, ecomb;
    Chi2Fit(8, sig, esig, comb, ecomb);
    if( TMath::Finite(comb+ecomb) )
      qt_all->Fill(ipt, 2*4*3+nstart, xpt, comb, ecomb);
  } //ipt

  mc(16+nstart, 2,2);
  legi(0, 0.2,0.8,0.7,0.9);
  leg0->SetNColumns(2);
  leg0->SetTextSize(0.02);

  for(int ibg=0; ibg<3; ibg++)
    for(int icr=0; icr<2; icr++)
      for(int pattern=0; pattern<4; pattern++)
      {
        mcd(16+nstart, ibg+1);
        int igr = icr + 2*pattern + 2*4*ibg;
        TGraphErrors *gr_all = qt_all->Graph(igr+nstart);

        gr_all->SetTitle( Form("#pi^{0} A_{LL} %s",region[ibg]) );
        aset(gr_all, "p_{T} [GeV]","A_{LL}", 0.,20., -0.2,0.4);
        style(gr_all, 1, 1+igr%8);
        gr_all->SetMarkerSize(0.);
        if(igr%8 == 0)
          gr_all->Draw("AP");
        else
          gr_all->Draw("P");
        if(ibg == 0)
          leg0->AddEntry(gr_all, Form("%s %s",pattern_list[pattern],crossing_list[icr]), "L");
      } // ibg, icr, pattern
  leg0->Draw();

  mcd(16+nstart, 4);
  gPad->SetGridy();
  TGraphErrors *gr_all = qt_all->Graph(2*4*3+nstart);
  gr_all->SetTitle( Form("#pi^{0} A_{LL} %s",region[3]) );
  aset(gr_all, "p_{T} [GeV]","A_{LL}", 0.,20., -0.01,0.03);
  style(gr_all, 1, 1);
  gr_all->SetMarkerSize(0.);
  gr_all->Draw("AP");

  qt_all->Write();
  for(int ibg=0; ibg<2; ibg++)
    for(int icr=0; icr<2; icr++)
      for(int pattern=0; pattern<4; pattern++)
      {
        int igr = icr + 2*pattern + 2*4*ibg + nstart;
        mcw( igr, Form("beam%d-isotype%d-pttype%d-bg%d-pattern%d-cross%d",
              beam, isotype, pttype, ibg, pattern, icr) );
      } // ibg, icr, pattern
  mcw( 16+nstart, Form("beam%d-isotype%d-pttype%d-combined",
        beam, isotype, pttype) );
}

void draw_IsoPionALL()
{
  QueryTree *qt_all = new QueryTree("data/IsoPionALL.root", "RECREATE");

  cout.precision(4);
  for(int beam=0; beam<3; beam++)
    for(int isotype=0; isotype<2; isotype++)
      for(int pttype=0; pttype<2; pttype++)
        draw_Each(qt_all, beam, isotype, pttype);

  qt_all->Close();
}
