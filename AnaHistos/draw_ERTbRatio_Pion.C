#include "FitMinv.h"

void draw_ERTbRatio_Pion()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-TAXI/PhotonHistos-total.root");

  // h[evtype][part]
  TH2 *h2_pion[3][3];

  int bbc10cm = 1;
  int tof = 1;
  int prob = 1;
  int checkmap = 1;
  int ival = 1;

  TH2 *h2_pion_t = (TH2*)f->Get("h2_pion_0");
  h2_pion_t = (TH2*)h2_pion_t->Clone();
  h2_pion_t->Reset();

  for(int evtype=0; evtype<3; evtype++)
    for(int part=0; part<3; part++)
    {
      h2_pion[evtype][part] = (TH2*)h2_pion_t->Clone(Form("h2_pion_type%d_part%d",evtype,part));
      for(int evenodd=0; evenodd<2; evenodd++)
        for(int pattern=0; pattern<3; pattern++)
          for(int isolated=0; isolated<2; isolated++)
          {
            int ih = part + 3*evenodd + 3*2*pattern + 3*2*3*evtype + 3*2*3*4*tof + 3*2*3*4*2*prob + 3*2*3*4*2*2*bbc10cm + 3*2*3*4*2*2*2*checkmap + 3*2*3*4*2*2*2*2*isolated + 3*2*3*4*2*2*2*2*2*ival;
            TH2 *h2_tmp = (TH2*)f->Get(Form("h2_pion_%d",ih));
            h2_pion[evtype][part]->Add(h2_tmp);
            delete h2_tmp;
          }
    }

  mc(0, 3,2);

  for(int part=0; part<3; part++)
  {
    double npion_ertc, enpion_ertc;
    double npion_ertb, enpion_ertb;
    TH1 *h_minv;

    mcd(0, part+1);
    h_minv = (TH1*)h2_pion[2][part]->ProjectionY("h_py", 23,30)->Clone("h_minv");
    h_minv->SetTitle( Form("ERT4x4c Part %d",part) );
    FitMinv(h_minv, npion_ertc, enpion_ertc, false, 0.10,0.17);
    delete h_minv;

    mcd(0, part+4);
    h_minv = (TH1*)h2_pion[1][part]->ProjectionY("h_py", 23,30)->Clone("h_minv");
    h_minv->SetTitle( Form("ERT4x4b Part %d",part) );
    FitMinv(h_minv, npion_ertb, enpion_ertb, false, 0.10,0.17);
    delete h_minv;

    double ratio = npion_ertc / npion_ertb;
    double eratio = ratio * sqrt( pow(enpion_ertc/npion_ertc,2) + pow(enpion_ertb/npion_ertb,2) );
    cout << "Part " << part << ", ratio = " << ratio << ", eratio = " << eratio << endl;
  }

  c0->Print("plots/ERTbRatio-pion.pdf");
}
