void AddHisto_Sasha(const int process = 0)
{
  const char *tname[4] = {"", "_t", "_p", "_tp"};

  TFile *f_out = new TFile(Form("histos/Pi0PP-%d.root",process), "RECREATE");

  TH1 *mchist[3][25][4];  // mchist[is][ipt][it]
  for(int is=0; is<3; is++)
    for(int ipt=0; ipt<25; ipt++)
      for(int it=0; it<4; it++)
      {
        char hname[100];
        sprintf(hname, "mchist_s%d_pt%02d%s", is, ipt, tname[it]);
        mchist[is][ipt][it] = new TH1F(hname, hname, 1000, 0., 1.);
      }

  int ptl[25];
  int ptr[25];
  int ipt = 0;
  for(int ip=0; ip<20; ip++)
  {
    ptl[ipt] = ip;
    ptr[ipt] = ip;
    ipt++;
  }
  for(int ip=20; ip<40; ip+=4)
  {
    ptl[ipt] = ip;
    ptr[ipt] = ip+3;
    ipt++;
  }

  const int nThread = 20;
  int thread = -1;
  int runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist-Sasha.txt");

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510ERT/12232/data/Pi0PP-%d.root",runnumber));
    //TFile *f = new TFile(Form("/phenix/spin/phnxsp01/shura/taxi/Run13pp510ERT/5094/data/%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *mch[3][40][4];  // mch[is][ip][it]
    for(int is=0; is<3; is++)
      for(int ip=0; ip<40; ip++)
        for(int it=0; it<4; it++)
        {
          char hname[100];
          sprintf(hname, "mc_s%d_bcc0_pt_%03d%s", is, 5*ip, tname[it]);
          mch[is][ip][it] = (TH1*)f->Get(hname);
        }

    for(int is=0; is<3; is++)
      for(int ipt=0; ipt<25; ipt++)
        for(int ip=ptl[ipt]; ip<=ptr[ipt]; ip++)
          for(int it=0; it<4; it++)
            mchist[is][ipt][it]->Add(mch[is][ip][it]);
  }

  f_out->cd();
  f_out->Write();
  f_out->Close();
}
