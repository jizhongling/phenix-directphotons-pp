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

    //TFile *f_mine = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    TFile *f_mine = new TFile(Form("/phenix/spin/phnxsp01/zji/taxi/Run13pp510ERT/13527/data/PhotonHistos-%d.root",runnumber));
    TFile *f_sasha = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510ERT/12232/data/Pi0PP-%d.root",runnumber));
    if( f_mine->IsZombie() || f_sasha->IsZombie() ) continue;

    for(int part=0; part<3; part++)
    {
      npion_mine[part] = 0.;
      npion_sasha[part] = 0.;
    }

    THnSparse *hn_pion = (THnSparse*)f_mine->Get("hn_pion");
    TAxis *axis_sec = hn_pion->GetAxis(0);
    TAxis *axis_pt = hn_pion->GetAxis(1);
    TAxis *axis_minv = hn_pion->GetAxis(2);
    TAxis *axis_pattern = hn_pion->GetAxis(3);
    TAxis *axis_cut = hn_pion->GetAxis(4);
    TAxis *axis_type = hn_pion->GetAxis(5);
    TAxis *axis_bbc10cm = hn_pion->GetAxis(6);

    for(int part=0; part<3; part++)
    {
      axis_bbc10cm->SetRange(2,2);
      axis_type->SetRange(3,3);
      axis_cut->SetRange(4,4);
      axis_sec->SetRange(secl[part],sech[part]);
      axis_pt->SetRange(5,25);
      TH1 *h_minv = hn_pion->Projection(2);
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
