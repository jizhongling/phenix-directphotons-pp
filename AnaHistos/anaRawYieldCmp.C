void anaRawYieldCmp(const Int_t process = 0)
{
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  const Int_t nThread = 20;
  Int_t thread = -1;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  TTree *t1 = new TTree("t1", "Raw yield");
  Int_t runnumber;
  Double_t npion_mine[3];
  Double_t npion_sasha[3];
  t1->Branch("runnumber", &runnumber, "runnumber/I");
  t1->Branch("npion_mine", npion_mine, "npion_mine[3]/D");
  t1->Branch("npion_sasha", npion_sasha, "npion_sasha[3]/D");

  TGraph *gr[3];
  Int_t igp[3] = {};
  for(Int_t part=0; part<3; part++)
  {
    gr[part] = new TGraph(nThread);
    gr[part]->SetName(Form("gr_%d",part));
  }

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f_mine = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    //TFile *f_sasha = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510ERT/12232/data/Pi0PP-%d.root",runnumber));
    TFile *f_sasha = new TFile(Form("/phenix/spin/phnxsp01/shura/taxi/Run13pp510ERT/5094/data/%d.root",runnumber));
    if( f_mine->IsZombie() || f_sasha->IsZombie() ) continue;

    for(Int_t part=0; part<3; part++)
    {
      npion_mine[part] = 0.;
      npion_sasha[part] = 0.;
    }

    THnSparse *hn_pion = (THnSparse*)f_mine->Get("hn_pion");
    for(Int_t part=0; part<3; part++)
    {
      hn_pion->GetAxis(3)->SetRange(2,2);
      hn_pion->GetAxis(0)->SetRange(secl[part],sech[part]);
      hn_pion->GetAxis(1)->SetRange(5,25);
      TH1 *h_minv = hn_pion->Projection(2);
      npion_mine[part] = h_minv->Integral(113,162);
      delete h_minv;
    }

    TH1 *mchist[3][40];  // mchist[is][ip]
    for(Int_t is=0; is<3; is++)
      for(Int_t ip=0; ip<40; ip++)
        mchist[is][ip] = (TH1*)f_sasha->Get(Form("mc_s%d_bcc0_pt_%03d_tp",is,5*ip));

    for(Int_t is=0; is<3; is++)
      for(Int_t ip=4; ip<40; ip++)
        npion_sasha[is] += mchist[is][ip]->Integral(113,162);

    t1->Fill();
    for(Int_t part=0; part<3; part++)
    {
      Double_t ratio = npion_mine[part] / npion_sasha[part];
      gr[part]->SetPoint(igp[part], runnumber, ratio);
      igp[part]++;
    }
    delete f_mine;
    delete f_sasha;
  }

  TFile *f_out = new TFile(Form("histos/RawYieldCmp-%d.root",process), "RECREATE");
  t1->Write();
  for(Int_t part=0; part<3; part++)
  {
    gr[part]->Set(igp[part]);
    gr[part]->Write();
  }
  f_out->Close();
}
