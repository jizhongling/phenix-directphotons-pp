#include "QueryTree.h"

void anaDCCheck(const int process = 0)
{
  const int nThread = 20;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  QueryTree *qt_dc = new QueryTree(Form("histos/DCCheck-%d.root",process), "RECREATE");

  TH2 *h2_phi[2];
  h2_phi[0] = new TH2F("h2_phi_zp", "DC phi vs run;runnumber;#phi [rad];", 900,0.,900., 50,-1.,4.);
  h2_phi[1] = (TH2*)h2_phi[0]->Clone("h2_phi_zm");

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/13912/data/PhotonHistos-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *h_events = (TH1*)f->Get("h_events");
    TH3 *h3_live = (TH3*)f->Get("h3_dclive_0");
    h3_live = (TH3*)h3_live->Clone();
    h3_live->Reset();
    for(int qual=4; qual<64; qual++)
    {
      int ih = qual;
      TH3 *h3_tmp = (TH3*)f->Get( Form("h3_dclive_%d",ih) );
      h3_live->Add(h3_tmp);
      delete h3_tmp;
    }

    double nev = h_events->GetBinContent( h_events->GetXaxis()->FindBin("ert_c_30cm") );

    for(int ns=0; ns<2; ns++)
      for(int we=0; we<2; we++)
      {
        int ig = we + 2*ns;
        double nyield = h3_live->Integral(101-100*ns,200-100*ns, 1+25*we,25+25*we, 0,-1);
        double yy = nyield / nev;
        double eyy = yy * sqrt( 1./nyield + 1./nev );
        qt_dc->Fill(runnumber, ig, (double)runnumber, yy, eyy);
      }

    for(int ns=0; ns<2; ns++)
    {
      h3_live->GetXaxis()->SetRange(101-100*ns,200-100*ns);
      TH1 *h_phi = h3_live->Project3D("y");
      for(int bin=1; bin<=h_phi->GetNbinsX(); bin++)
        h2_phi[ns]->SetBinContent( thread+1, bin, h_phi->GetBinContent(bin)/nev );
    }

    delete h3_live;
    delete f;
  }

  qt_dc->Write();
  for(int ns=0; ns<2; ns++)
    h2_phi[ns]->Write();
  qt_dc->Close();
}
