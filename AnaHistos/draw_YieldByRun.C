void draw_YieldByRun()
{
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  Int_t runnumber;
  Double_t npion_notof[3], npion_tof[3];

  TTree *t_ert = new TTree("t_ert", "Pi0 raw yield from ERT sample");
  t_ert->Branch("runnumber", &runnumber, "runnumber/I");
  t_ert->Branch("npion_notof", npion_notof, "npion_notof[3]/D");
  t_ert->Branch("npion_tof", npion_tof, "npion_tof[3]/D");

  TTree *t_mb = new TTree("t_mb", "Pi0 raw yield from MB sample");
  t_mb->Branch("runnumber", &runnumber, "runnumber/I");
  t_mb->Branch("npion_notof", npion_notof, "npion_notof[3]/D");
  t_mb->Branch("npion_tof", npion_tof, "npion_tof[3]/D");

  ifstream fin_ert("/phenix/plhf/zji/taxi/Run13pp510ERT/runnumber.txt");
  while( fin_ert >> runnumber )
  {
    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ERT/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<3; is++)
      {
        hn_pion->GetAxis(3)->SetRange(ic+1,ic+1);
        hn_pion->GetAxis(0)->SetRange(secl[is],sech[is]);
        hn_pion->GetAxis(1)->SetRange(5,30);
        TH1 *h_minv = hn_pion->Projection(2);

        if(ic==0)
          npion_notof[is] = h_minv->Integral(113,162);
        else if(ic==1)
          npion_tof[is] = h_minv->Integral(113,162);

        delete h_minv;
      }

    t_ert->Fill();
    delete f;
  }
  fin_ert.close();

  ifstream fin_mb("/phenix/plhf/zji/taxi/Run13pp510MinBias/runnumber.txt");
  while( fin_mb >> runnumber )
  {
    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-MB/PhotonNode-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    THnSparse *hn_pion = (THnSparse*)f->Get("hn_pion");

    for(Int_t ic=0; ic<2; ic++)
      for(Int_t is=0; is<3; is++)
      {
        hn_pion->GetAxis(3)->SetRange(ic+1,ic+1);
        hn_pion->GetAxis(0)->SetRange(secl[is],sech[is]);
        hn_pion->GetAxis(1)->SetRange(5,30);
        TH1 *h_minv = hn_pion->Projection(2);

        if(ic==0)
          npion_notof[is] = h_minv->Integral(113,162);
        else if(ic==1)
          npion_tof[is] = h_minv->Integral(113,162);

        delete h_minv;
      }

    t_mb->Fill();
    delete f;
  }
  fin_mb.close();

  TFile *f_out = new TFile("data/RawYield.root", "RECREATE");
  t_ert->Write();
  t_mb->Write();
  f_out->Close();
}
