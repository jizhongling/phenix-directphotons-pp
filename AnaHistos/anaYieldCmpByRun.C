void anaYieldCmpByRun(const int process = 0)
{
  const int secl[3] = {1, 5, 7};
  const int sech[3] = {4, 6, 8};

  const int nThread = 20;
  int thread = -1;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runnumber.txt");

  TTree *t1 = new TTree("t1", "Raw yield");
  int runnumber;
  double npion_mine[3];
  double npion_sasha[3];
  t1->Branch("runnumber", &runnumber, "runnumber/I");
  t1->Branch("npion_mine", npion_mine, "npion_mine[3]/D");
  t1->Branch("npion_sasha", npion_sasha, "npion_sasha[3]/D");

  TGraph *gr[3];
  int igp[3] = {};
  for(int part=0; part<3; part++)
  {
    gr[part] = new TGraph(nThread);
    gr[part]->SetName(Form("gr_%d",part));
  }

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f_mine = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/13597/data/PhotonHistos-%d.root",runnumber));
    TFile *f_sasha = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510ERT/12232/data/Pi0PP-%d.root",runnumber));
    if( f_mine->IsZombie() || f_sasha->IsZombie() ) continue;

    for(int part=0; part<3; part++)
    {
      npion_mine[part] = 0.;
      npion_sasha[part] = 0.;
    }

    // h2_pion[part]
    TH2 *h2_pion[3];
    TH2 *h2_pion_t = (TH2*)f_mine->Get("h2_pion_0");
    h2_pion_t->Reset();
    int bbc10cm = 1;
    int evtype = 2;
    int cut = 3;
    for(int part=0; part<3; part++)
    {
      h2_pion[part] = (TH2*)h2_pion_t->Clone(Form("h2_pion_%d",part));
      for(int sector=secl[part]-1; sector<=sech[part]-1; sector++)
        for(int pattern=0; pattern<3; pattern++)
        {
          int ih = sector + 8*pattern + 3*8*cut + 4*3*8*evtype + 3*4*3*8*bbc10cm;
          TH2 *h2_tmp = (TH2*)f_mine->Get(Form("h2_pion_%d",ih));
          h2_pion[part]->Add(h2_tmp);
          delete h2_tmp;
        }
    }

    for(int part=0; part<3; part++)
    {
      TH1 *h_minv = h2_pion[part]->ProjectionY("h_minv", 5,25);
      npion_mine[part] = h_minv->Integral(113,162);
      delete h_minv;
    }

    TH1 *mchist[3][40];  // mchist[is][ip]
    for(int is=0; is<3; is++)
      for(int ip=0; ip<40; ip++)
        mchist[is][ip] = (TH1*)f_sasha->Get(Form("mc_s%d_bcc0_pt_%03d_tp",is,5*ip));

    for(int is=0; is<3; is++)
      for(int ip=4; ip<40; ip++)
        npion_sasha[is] += mchist[is][ip]->Integral(113,162);

    t1->Fill();
    for(int part=0; part<3; part++)
    {
      double ratio = npion_mine[part] / npion_sasha[part];
      gr[part]->SetPoint(igp[part], runnumber, ratio);
      igp[part]++;
    }
    delete f_mine;
    delete f_sasha;
  }

  TFile *f_out = new TFile(Form("histos/YieldCmpByRun-%d.root",process), "RECREATE");
  t1->Write();
  for(int part=0; part<3; part++)
  {
    gr[part]->Set(igp[part]);
    gr[part]->Write();
  }
  f_out->Close();
}
