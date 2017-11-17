void anaProbEff(const Int_t process = 0)
{
  TFile *f_out = new TFile(Form("eff/ProbEff-%d.root",process), "RECREATE");
  TH1 *h_pass[2][25];  // h_pass[part][ipt]
  TH1 *h_total[2][25];  // h_total[part][ipt]
  for(Int_t part=0; part<2; part++)
    for(Int_t ipt=0; ipt<25; ipt++)
    {
      char hname[100];
      sprintf(hname,"h_pass_part%d_pt%d",part,ipt);
      h_pass[part][ipt] = new TH1F(hname,hname,1000,0.,1.);
      sprintf(hname,"h_total_part%d_pt%d",part,ipt);
      h_total[part][ipt] = new TH1F(hname,hname,1000,0.,1.);
    }

  Int_t ptl[25];
  Int_t ptr[25];
  Int_t ipt = 0;
  for(Int_t ip=0; ip<20; ip++)
  {
    ptl[ipt] = ip;
    ptr[ipt] = ip;
    ipt++;
  }
  for(Int_t ip=20; ip<40; ip+=4)
  {
    ptl[ipt] = ip;
    ptr[ipt] = ip+3;
    ipt++;
  }

  const Int_t nThread = 20;
  Int_t thread = -1;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510ERT/12232/data/Pi0PP-%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *mc_pass[3][40];  // mc_pass[is][ip]
    TH1 *mc_total[3][40];  // mc_total[is][ip]
    for(Int_t is=0; is<3; is++)
      for(Int_t ip=0; ip<40; ip++)
      {
        char hname[100];
        sprintf(hname,"mc_s%d_bcc0_pt_%03d_tp",is,5*ip);
        mc_pass[is][ip] = (TH1*)f->Get(hname);
        sprintf(hname,"mc_s%d_bcc0_pt_%03d_t",is,5*ip);
        mc_total[is][ip] = (TH1*)f->Get(hname);
      }

    for(Int_t is=0; is<3; is++)
      for(Int_t ipt=0; ipt<25; ipt++)
        for(Int_t ip=ptl[ipt]; ip<=ptr[ipt]; ip++)
        {
          h_pass[is/2][ipt]->Add(mc_pass[is][ip]);
          h_total[is/2][ipt]->Add(mc_total[is][ip]);
        }
  }

  f_out->cd();
  f_out->Write();
  f_out->Close();
}
