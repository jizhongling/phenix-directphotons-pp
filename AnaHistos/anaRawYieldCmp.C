void anaRawYieldCmp(const Int_t process = 0)
{
  const Int_t secl[3] = {1, 5, 7};
  const Int_t sech[3] = {4, 6, 8};

  TGraph *gr[3];
  for(Int_t part=0; part<3; part++)
  {
    gr[part] = new TGraph(20);
    gr[part]->SetName(Form("gr_%d",part));
  }

  const Int_t nThread = 20;
  Int_t thread = -1;
  Int_t irun = 0;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f_mine = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/PhotonNode-%d.root",runnumber));
    TFile *f_sasha = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510MinBias/11955/data/Pi0PP-%d.root",runnumber));
    if( f_mine->IsZombie() || f_sasha->IsZombie() ) continue;

    Double_t npion_mine[3];
    Double_t npion_sasha[3];
    for(Int_t part=0; part<3; part++)
    {
      npion_mine[part] = 0.;
      npion_sasha[part] = 0.;
    }

    THnSparse *hn_pion = (THnSparse*)f_mine->Get("hn_pion");
    hn_pion->GetAxis(1)->SetRange(5,25);
    TH1 *h_minv[3];
    for(Int_t part=0; part<3; part++)
    {
      hn_pion->GetAxis(0)->SetRange(secl[part],sech[part]);
      h_minv[part] = hn_pion->Projection(2);
      h_minv[part]->SetName(Form("h_minv_%d",part));
    }

    Int_t bin047 = hn_pion->GetAxis(2)->FindBin(0.047);
    Int_t bin097 = hn_pion->GetAxis(2)->FindBin(0.097);
    Int_t bin112 = hn_pion->GetAxis(2)->FindBin(0.112);
    Int_t bin162 = hn_pion->GetAxis(2)->FindBin(0.162);
    Int_t bin187 = hn_pion->GetAxis(2)->FindBin(0.187);
    Int_t bin227 = hn_pion->GetAxis(2)->FindBin(0.227);

    for(Int_t part=0; part<3; part++)
      npion_mine[part] += h_minv[part]->Integral(bin112,bin162) - ( h_minv[part]->Integral(bin047,bin097) + h_minv[part]->Integral(bin187,bin227) ) / 2.;

    TH1 *mchist[3][40];  // mchist[is][ip]
    for(Int_t is=0; is<3; is++)
      for(Int_t ip=0; ip<40; ip++)
        mchist[is][ip] = (TH1*)f_sasha->Get(Form("mc_s%d_bcc0_pt_%03d_tp",is,5*ip));

    bin047 = mchist[0][0]->GetXaxis()->FindBin(0.047);
    bin097 = mchist[0][0]->GetXaxis()->FindBin(0.097);
    bin112 = mchist[0][0]->GetXaxis()->FindBin(0.112);
    bin162 = mchist[0][0]->GetXaxis()->FindBin(0.162);
    bin187 = mchist[0][0]->GetXaxis()->FindBin(0.187);
    bin227 = mchist[0][0]->GetXaxis()->FindBin(0.227);

    for(Int_t is=0; is<3; is++)
      for(Int_t ip=4; ip<40; ip++)
        npion_sasha[is] += mchist[is][ip]->Integral(bin112,bin162) - ( mchist[is][ip]->Integral(bin047,bin097) - mchist[is][ip]->Integral(bin187,bin227) ) / 2.;

    for(Int_t part=0; part<3; part++)
    {
      Double_t ratio = npion_mine[part] / npion_sasha[part];
      gr[part]->SetPoint(irun, runnumber, ratio);
      delete h_minv[part];
    }
    delete f_mine;
    delete f_sasha;
    irun++;
  }

  TFile *f_out = new TFile(Form("histos/RawYieldCmp-%d.root",process), "RECREATE");
  for(Int_t part=0; part<3; part++)
    gr[part]->Write();
  f_out->Close();
}
