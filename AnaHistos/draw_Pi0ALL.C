#include "GlobalVars.h"
#include "QueryTree.h"
#include "Ftest.h"

void draw_Pi0ALL()
{
  const char *crossing_list[2] = {"even", "odd"};
  const char *pattern_list[4] = {"SOOSSOO", "OSSOOSS", "SSOO", "OOSS"};
  const char *region[2] = {"Sig+BG", "BG"};

  QueryTree *qt_all = new QueryTree("data/Pi0ALL.root", "RECREATE");

  QueryTree *qt_asym = new QueryTree("data/pion-asym.root");
  qt_asym->SetQuiet();

  vector<double> *vp_ALL = new vector<double>[8];
  vector<double> *vp_eALL = new vector<double>[8];

  for(int ibg=0; ibg<2; ibg++)
    for(int ipt=0; ipt<npT_pol; ipt++)
    {
      for(int icr=0; icr<2; icr++)
        for(int pattern=0; pattern<4; pattern++)
        {
          int id = icr + 2*pattern;
          vp_ALL[id].clear();
          vp_eALL[id].clear();

          int ig = ibg + 2*icr + 2*2*pattern + 2*2*4*ipt;
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
      cout << pTbin_pol[ipt] << "-" << pTbin_pol[ipt+1] << ": " << F << ", " << p << endl;
    } // ibg, ipt

  TF1 *fn_mean = new TF1("fn_mean", "pol0");

  for(int ibg=0; ibg<2; ibg++)
    for(int icr=0; icr<2; icr++)
      for(int pattern=0; pattern<4; pattern++)
      {
        int igr = icr + 2*pattern + 2*4*ibg;
        mc(igr);

        for(int ipt=0; ipt<npT_pol; ipt++)
        {
          int ig = ibg + 2*icr + 2*2*pattern + 2*2*4*ipt;

          mcd(igr, ipt+1);
          TGraphErrors *gr = qt_asym->Graph(ig); 
          gr->SetTitle(Form("p_{T}: %.1f-%.1f GeV",pTbin_pol[ipt],pTbin_pol[ipt+1]));
          aset(gr, "runnumber","A_{LL}", 386700.,398200., -1.,1.);
          style(gr, 20, 1);
          gr->Draw("AP");  // must before GetXaxis()

          gr->Fit(fn_mean, "Q");
          double mean = fn_mean->GetParameter(0);
          double emean = fn_mean->GetParError(0);

          if( TMath::Finite(mean+emean) )
          {
            double xpt = ( pTbin_pol[ipt] + pTbin_pol[ipt+1] ) / 2.;
            qt_all->Fill(ipt, igr, xpt, mean, emean);
          }
        } // ipt
      } // ibg, icr, pattern

  mc(16, 2);
  legi(0, 0.2,0.8,0.7,0.9);
  leg0->SetNColumns(2);
  leg0->SetTextSize(0.02);
  for(int ibg=0; ibg<2; ibg++)
    for(int icr=0; icr<2; icr++)
      for(int pattern=0; pattern<4; pattern++)
      {
        mcd(16, ibg+1);
        int igr = icr + 2*pattern + 2*4*ibg;
        TGraphErrors *gr_all = qt_all->Graph(igr);

        gr_all->SetTitle( Form("#pi^{0} A_{LL} %s",region[ibg]) );
        aset(gr_all, "p_{T} [GeV]","A_{LL}", 0.,30., -1.,1.);
        style(gr_all, 1, 1+igr%8);
        gr_all->SetMarkerSize(0.);
        if(igr%8==0)
          gr_all->Draw("AP");
        else
          gr_all->Draw("P");
        if(ibg==0)
          leg0->AddEntry(gr_all, Form("%s %s",pattern_list[pattern],crossing_list[icr]), "L");
      } // ibg, icr, pattern
  leg0->Draw();
  c16->Print("plots/Pi0ALL.pdf");

  qt_all->Write();
  for(int ibg=0; ibg<2; ibg++)
    for(int icr=0; icr<2; icr++)
      for(int pattern=0; pattern<4; pattern++)
      {
        int igr = icr + 2*pattern + 2*4*ibg;
        mcw( igr, Form("bg%d-pattern%d-cross%d", ibg, pattern, icr) );
      } // ibg, icr, pattern
  qt_all->Close();
}
