#include "QueryTree.h"
#include "MultiGraph.h"

void draw_DCCheck(const int print_dcboard = 0)
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h2/h3[ns][we]
  TH3 *h3_dphiz[2][2];
  TH2 *h2_board[2][2];

  TH3 *h3_dphiz_t = (TH3*)f->Get("h3_dcdphiz_0");
  TH2 *h2_board_t = (TH2*)f->Get("h2_alphaboard_0");
  TH3 *h3_live = (TH3*)f->Get("h3_dclive_0");
  h3_dphiz_t = (TH3*)h3_dphiz_t->Clone();
  h2_board_t = (TH2*)h2_board_t->Clone();
  h3_live = (TH3*)h3_live->Clone();
  h3_dphiz_t->Reset();
  h2_board_t->Reset();
  h3_live->Reset();

  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      h3_dphiz[ns][we] = (TH3*)h3_dphiz_t->Clone( Form("h3_dphiz_ns%d_we_%d",ns,we) );
      h2_board[ns][we] = (TH2*)h2_board_t->Clone( Form("h2_board_ns%d_we_%d",ns,we) );
      for(int iqual=1; iqual<3; iqual++)
      {
        int ih = ns + 2*we + 2*2*iqual;
        TH3 *h3_dphiz_tmp = (TH3*)f->Get( Form("h3_dcdphiz_%d",ih) );
        TH2 *h2_board_tmp = (TH2*)f->Get( Form("h2_alphaboard_%d",ih) );
        h3_dphiz[ns][we]->Add(h3_dphiz_tmp);
        h2_board[ns][we]->Add(h2_board_tmp);
        delete h3_dphiz_tmp;
        delete h2_board_tmp;
      }
    }

  for(int iqual=1; iqual<3; iqual++)
  {
    int ih = iqual;
    TH3 *h3_tmp = (TH3*)f->Get( Form("h3_dclive_%d",ih) );
    h3_live->Add(h3_tmp);
    delete h3_tmp;
  }

  mc(0, 2,2);
  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      mcd(0, 2*ns+we+1);
      TString NS = ns ? "S" : "N";
      TString WE = we ? "E" : "W";
      h2_board[ns][we]->SetTitle(NS+WE);
      h2_board[ns][we]->DrawCopy("COLZ");
    }

  mc(1);
  mcd(1);
  h3_live->GetZaxis()->SetRange(3,-1);
  h3_live->Project3D("yx")->Draw("COLZ");

  TF1 *f_fit = new TF1("f_fit", "gaus(0)+pol2(3)", -100., 100.);
  for(int id=0; id<2; id++)
  {
    mc(id+2, 2,2);
    for(int ns=0; ns<2; ns++)
      for(int we=0; we<2; we++)
      {
        TString ax = id ? "y" : "x";
        TString NS = ns ? "S" : "N";
        TString WE = we ? "E" : "W";

        mcd(id+2, 2*ns+we+1);
        h3_dphiz[ns][we]->GetZaxis()->SetRange(5,15);
        TH1 *h_diff = h3_dphiz[ns][we]->Project3D(ax);
        double max = h_diff->GetMaximum();
        double sigma = id ? 3.3 : 0.005;
        double par[10] = {max,0.,sigma, 0.,0.,0.};
        f_fit->SetParameters(par);
        h_diff->Fit(f_fit, "Q0","", -5.*sigma,5.*sigma);
        sigma = f_fit->GetParameter(2);
        h_diff->GetXaxis()->SetRangeUser(-3.*sigma,3.*sigma);
        cout << NS+WE+" RMS: " << h_diff->GetRMS() << endl;

        h_diff->GetXaxis()->SetRange(0,-1);
        h_diff->SetTitle(NS+WE);
        aset(h_diff);
        f_fit->SetLineColor(kRed);
        h_diff->DrawCopy();
        f_fit->DrawCopy("SAME");
      }
  }

  QueryTree *qt_dc = new QueryTree("data/DCCheck.root");

  mc(4, 2,2);
  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      int ig = we + 2*ns;
      mcd(4, ig+1);
      TString NS = ns ? "S" : "N";
      TString WE = we ? "E" : "W";
      TGraphErrors *gr = qt_dc->Graph(ig);
      gr->SetTitle(NS+WE);
      aset(gr, "runnumber","Ntrack/Nevent", 387000.,398200., 0.,1.5);
      style(gr, 20, 1);
      gr->Draw("AP");

      double mean, sigma;
      GetMeanSigma<TGraphErrors>(gr, mean, sigma);
      double ylow = mean - 3. * sigma;
      TLine *line = new TLine;
      line->SetLineColor(kRed);
      line->DrawLine(387000.,mean,398200.,mean);
      line->SetLineColor(kGreen);
      line->DrawLine(387000.,ylow,398200.,ylow);

      int N = gr->GetN();
      for(int i=0; i<N; i++)
      {
        double xx, yy;
        gr->GetPoint(i, xx, yy);
        if( yy < ylow )
          cout << xx << " ";
      }
    }
  cout << endl;

  mc(5, 2,1);
  TH2 *h2_phi[2];  // h2_phi[ns]
  for(int ns=0; ns<2; ns++)
  {
    mcd(5, ns+1);
    TString name = ns ? "h2_phi_zm" : "h2_phi_zp";
    TString title = ns ? "DC phi vs run in South" : "DC phi vs run in North";
    h2_phi[ns] = (TH2*)qt_dc->Get(name);
    h2_phi[ns]->SetTitle(title);
    h2_phi[ns]->Draw("COLZ");
  }

  if(!print_dcboard)
    return;

  mc(6, 2,2);
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist-DC3sigma.txt");
  c6->Print("plots/DCAlphaBoard.pdf(", "pdf");
  while( fin >> runnumber )
  {
    for(int ns=0; ns<2; ns++)
      for(int we=0; we<2; we++)
      {
        mcd(6, ns*2+we+1);
        TString NS = ns ? "S" : "N";
        TString WE = we ? "E" : "W";
        TH2 *h2_board_run = (TH2*)qt_dc->Get( Form("h2_board_ns%d_we%d_%d",ns,we,runnumber) );
        if(!h2_board_run) continue;
        h2_board_run->SetTitle( Form("%d, %s",runnumber,(NS+WE).Data()) );
        h2_board_run->DrawCopy("COLZ");
        delete h2_board_run;
      }
    c6->Print("plots/DCAlphaBoard.pdf", "pdf");
    c6->Clear("D");
  }
  c6->Print("plots/DCAlphaBoard.pdf)", "pdf");
}
