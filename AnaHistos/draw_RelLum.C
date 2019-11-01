void draw_RelLum()
{
  const char *crname[2] = {"even", "odd"};
  const char *rname[4] = {"GL1p", "SS. Uncorr", "SS. Pileup", "SS. Residual"};

  TFile *f_rlum = new TFile("data/RelLum.root");
  TTree *t_rlum_check = (TTree*)f_rlum->Get("T1");
  int runnumber;
  bool runqa;
  double rlum[4][2], erlum[4][2];
  t_rlum_check->SetBranchAddress("Runnumber", &runnumber);
  t_rlum_check->SetBranchAddress("RunQA", &runqa);
  t_rlum_check->SetBranchAddress("GL1p", rlum[0]);
  t_rlum_check->SetBranchAddress("eGL1p", erlum[0]);
  t_rlum_check->SetBranchAddress("SS_Uncorr", rlum[1]);
  t_rlum_check->SetBranchAddress("eSS_Uncorr", erlum[1]);
  t_rlum_check->SetBranchAddress("SS_Pileup", rlum[2]);
  t_rlum_check->SetBranchAddress("eSS_Pileup", erlum[2]);
  t_rlum_check->SetBranchAddress("SS_Residual", rlum[3]);
  t_rlum_check->SetBranchAddress("eSS_Residual", erlum[3]);

  int nentries = t_rlum_check->GetEntries();

  TGraphErrors *gr_rlum[4][2], *gr_ratio[3][2];
  int igr = 0;
  for(int ir=0; ir<4; ir++)
    for(int icr=0; icr<2; icr++)
    {
      gr_rlum[ir][icr] = new TGraphErrors(nentries);
      if(ir > 0)
        gr_ratio[ir-1][icr] = new TGraphErrors(nentries);
    }

  for(int ien=0; ien<nentries; ien++)
  {
    t_rlum_check->GetEntry(ien);
    if(!runqa) continue;

    for(int ir=0; ir<4; ir++)
      for(int icr=0; icr<2; icr++)
      {
        gr_rlum[ir][icr]->SetPoint(igr, (double)runnumber, rlum[ir][icr]);
        gr_rlum[ir][icr]->SetPointError(igr, 0., erlum[ir][icr]);
        if(ir > 0)
        {
          double ratio = rlum[ir][icr]/rlum[0][icr];
          double eratio = ratio*sqrt(pow(erlum[ir][icr]/rlum[ir][icr],2) + pow(erlum[0][icr]/rlum[0][icr],2));
          gr_ratio[ir-1][icr]->SetPoint(igr, (double)runnumber, ratio);
          gr_ratio[ir-1][icr]->SetPointError(igr, 0., eratio);
        }
      } // ir, icr
    igr++;
  } // ien

  mc(0, 2,2);
  legi(0, 0.6,0.7,0.9,0.9);
  legi(1, 0.2,0.7,0.5,0.9);

  for(int icr=0; icr<2; icr++)
  {
    mcd(0, icr+1);
    for(int ir=0; ir<4; ir++)
    {
      gr_rlum[ir][icr]->Set(igr);
      gr_rlum[ir][icr]->SetTitle( Form("Relative luminosity %s crossing", crname[icr]) );
      aset(gr_rlum[ir][icr], "Runnumber","Relative luminosity", 365000.,369000., 0.8,1.2);
      style(gr_rlum[ir][icr], 24+ir, 1+ir);
      if(ir == 0)
        gr_rlum[ir][icr]->Draw("AP");
      else 
        gr_rlum[ir][icr]->Draw("P");
      if(icr == 0)
        leg0->AddEntry(gr_rlum[ir][icr], rname[ir], "P");
    }
    leg0->Draw();

    mcd(0, icr+3);
    for(int ir=0; ir<3; ir++)
    {
      gr_ratio[ir][icr]->Set(igr);
      gr_ratio[ir][icr]->SetTitle( Form("Ratio of RelLum %s crossing", crname[icr]) );
      aset(gr_ratio[ir][icr], "Runnumber","Ratio of RelLum", 365000.,369000., 0.98,1.02);
      style(gr_ratio[ir][icr], 25+ir, 2+ir);
      if(ir==0)
        gr_ratio[ir][icr]->Draw("AP");
      else 
        gr_ratio[ir][icr]->Draw("P");
      if(icr == 0)
        leg1->AddEntry(gr_ratio[ir][icr], rname[ir+1], "P");
    }
    leg1->Draw();
  }

  c0->Print("plots/RelLum.pdf");

}
