#include "QueryTree.h"

void anaDCCheck(const int process = 0)
{
  gSystem->Load("libDirectPhotonPP.so");

  const int nThread = 20;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist-Sasha.txt");

  QueryTree *qt_dc = new QueryTree(Form("histos/DCCheck-%d.root",process), "RECREATE");

  TH2 *h2_phi[2];
  h2_phi[0] = new TH2F("h2_phi_zp", "DC phi vs run;runnumber;#phi [rad];", 900,0.,900., 50,-1.,4.);
  h2_phi[1] = (TH2*)h2_phi[0]->Clone("h2_phi_zm");

  DCDeadmapChecker *dcdeadmap = new DCDeadmapChecker();

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread )
      continue;
    else if( thread >= (process+1)*nThread )
      break;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/15811/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("h_events");
    TH2 *h2_board_t = (TH2*)f->Get("h2_alphaboard_0");
    TH3 *h3_live_t = (TH3*)f->Get("h3_dclive_0");
    h2_board_t = (TH2*)h2_board_t->Clone();
    h3_live_t = (TH3*)h3_live_t->Clone();
    h2_board_t->Reset();
    h3_live_t->Reset();

    dcdeadmap->SetMapByRunnumber(runnumber);

    TH2 *h2_board[2][2];
    for(int ns=0; ns<2; ns++)
      for(int we=0; we<2; we++)
      {
        h2_board[ns][we] = (TH2*)h2_board_t->Clone( Form("h2_board_ns%d_we%d_%d",ns,we,runnumber) );
        for(int iqual=1; iqual<3; iqual++)
        {
          int ih = ns + 2*we + 2*2*iqual;
          TH2 *h2_board_tmp = (TH2*)f->Get( Form("h2_alphaboard_%d",ih) );
          h2_board[ns][we]->Add(h2_board_tmp);
          delete h2_board_tmp;
        }
        string dcns = ns ? "S" : "N";
        string dcwe = we ? "E" : "W";
        string nswe = dcns + dcwe;
        for(int binx=1; binx<=h2_board[ns][we]->GetXaxis()->GetLast(); binx++)
          for(int biny=1; biny<=h2_board[ns][we]->GetYaxis()->GetLast(); biny++)
          {
            double board = h2_board[ns][we]->GetXaxis()->GetBinCenter(binx);
            double alpha = h2_board[ns][we]->GetYaxis()->GetBinCenter(biny);
            if( dcdeadmap->IsDead(nswe, board, alpha) )
              h2_board[ns][we]->SetBinContent(binx, biny, 0.);
          }
        qt_dc->cd();
        h2_board[ns][we]->Write();
      }

    TH3 *h3_live = (TH3*)h3_live_t->Clone("h3_live");
    for(int iqual=1; iqual<3; iqual++)
    {
      int ih = iqual;
      TH3 *h3_live_tmp = (TH3*)f->Get( Form("h3_dclive_%d",ih) );
      h3_live->Add(h3_live_tmp);
      delete h3_live_tmp;
    }

    double nev = h_events->GetBinContent( h_events->GetXaxis()->FindBin("ert_c_30cm") );

    for(int ns=0; ns<2; ns++)
      for(int we=0; we<2; we++)
      {
        int ig = we + 2*ns;
        //double nyield = h3_live->Integral(101-100*ns,200-100*ns, 1+25*we,25+25*we, 0,-1);
        double nyield = h2_board[ns][we]->Integral(0,-1, 0,-1);
        double yy = nyield / nev;
        double eyy = yy * sqrt( 1./nyield + 1./nev );
        if( TMath::Finite(yy+eyy) )
          qt_dc->Fill(runnumber, ig, (double)runnumber, yy, eyy);
        else
          qt_dc->Fill(runnumber, ig, (double)runnumber, 0., 1.);
      }

    for(int ns=0; ns<2; ns++)
    {
      h3_live->GetXaxis()->SetRange(101-100*ns,200-100*ns);
      TH1 *h_phi = h3_live->Project3D("y");
      for(int bin=1; bin<=h_phi->GetXaxis()->GetLast(); bin++)
        h2_phi[ns]->SetBinContent( thread+1, bin, h_phi->GetBinContent(bin)/nev );
    }

    delete h2_board_t;
    delete h3_live_t;
    delete h3_live;
    for(int ns=0; ns<2; ns++)
      for(int we=0; we<2; we++)
        delete h2_board[ns][we];
    delete f;
  }

  qt_dc->Write();
  for(int ns=0; ns<2; ns++)
    h2_phi[ns]->Write();
  qt_dc->Close();
  delete dcdeadmap;
}
