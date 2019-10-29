void anaSSkFactor()
{
  TFile *f_bbc = new TFile("data/BBC_Run13_PING.root");
  TTree *t_bbc = (TTree*)f_bbc->Get("T");
  int runno, clock[120], bbcn[120], bbcs[120], bbcns[120], bbc30[120], gl1p[120];
  t_bbc->SetBranchAddress("Runnumber", &runno);
  t_bbc->SetBranchAddress("CLOCK", clock);
  t_bbc->SetBranchAddress("BBCN", bbcn);
  t_bbc->SetBranchAddress("BBCS", bbcs);
  t_bbc->SetBranchAddress("BBCNS", bbcns);
  t_bbc->SetBranchAddress("BBC_30cm", bbc30);
  t_bbc->SetBranchAddress("GL1p_BBC30", gl1p);

  TFile *f_k = new TFile("data/SSkFactor.root", "RECREATE");
  TTree *t_k = new TTree("T", "k factor");
  double rns, rns_true, kn, ks, ekn, eks;
  t_k->Branch("Rate", &rns, "Rate/D");
  t_k->Branch("Rate_true", &rns_true, "Rate_true/D");
  t_k->Branch("kN", &kn, "kN/D");
  t_k->Branch("kS", &ks, "kS/D");
  t_k->Branch("ekN", &ekn, "ekN/D");
  t_k->Branch("ekS", &eks, "ekS/D");

  TGraphErrors *gr_qa = new TGraphErrors(120);

  for(int ien=0; ien<t_bbc->GetEntries(); ien++)
  {
    if(ien && ien%100 == 0)
      cout << ien << endl;

    t_bbc->GetEntry(ien);

    int igr = 0;
    for(int ib=0; ib<120; ib++)
      if( bbc30[ib] > 0 && gl1p[ib] > 0 )
      {
        double rbbc = (double)gl1p[ib]/bbc30[ib];
        double erbbc = rbbc*sqrt(1./gl1p[ib] + 1./bbc30[ib]);
        gr_qa->SetPoint(igr, (double)ib, rbbc);
        gr_qa->SetPointError(igr, 0., erbbc);
        igr++;
      }
    gr_qa->Set(igr);

    if(igr == 0) continue;
    TFitResultPtr r_qa = gr_qa->Fit("pol0", "QS");
    double par0 = r_qa->Value(0);
    double chi2 = r_qa->Chi2();
    unsigned ndf = r_qa->Ndf();
    if( par0 < 0.998 || par0 > 1.002 ||
        ndf == 0 || chi2/ndf > 2.5e3 )
      continue;

    //cout << runno << " ";

    for(int ib=0; ib<120; ib++)
      if( clock[ib] > 0 && bbcn[ib] > 0 && bbcs[ib] > 0 && bbcns[ib] > 0 )
      {
        rns = (double)bbcns[ib]/clock[ib];
        double rn = (double)bbcn[ib]/clock[ib];
        double rs = (double)bbcs[ib]/clock[ib];

        double erns2 = rns*rns*(1./bbcns[ib] + 1./clock[ib]);
        double ern2 = rn*rn*(1./bbcn[ib] + 1./clock[ib]);
        double ers2 = rs*rs*(1./bbcs[ib] + 1./clock[ib]);

        rns_true = log(1 - rn - rs + rns) - log(1 - rn) - log(1 - rs);
        double rn_true = -log(1 - rn) - rns_true;
        double rs_true = -log(1 - rs) - rns_true;

        kn = rn_true/rns_true;
        ks = rs_true/rns_true;

        ekn = sqrt(ers2*pow(-(((1/(1 - rs) - 1/(1 - rn + rns - rs))*(log(1 - rs) - log(1 - rn + rns - rs)))/pow(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs),2)) + (-(1/(1 - rs)) + 1/(1 - rn + rns - rs))/(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs)),2) + erns2*pow(-((log(1 - rs) - log(1 - rn + rns - rs))/((1 - rn + rns - rs)*pow(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs),2))) - 1/((1 - rn + rns - rs)*(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs))),2) + ern2*pow(-(((1/(1 - rn) - 1/(1 - rn + rns - rs))*(log(1 - rs) - log(1 - rn + rns - rs)))/pow(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs),2)) + 1/((1 - rn + rns - rs)*(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs))),2));
        eks = sqrt(ern2*pow(-(((1/(1 - rn) - 1/(1 - rn + rns - rs))*(log(1 - rn) - log(1 - rn + rns - rs)))/pow(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs),2)) + (-(1/(1 - rn)) + 1/(1 - rn + rns - rs))/(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs)),2) + erns2*pow(-((log(1 - rn) - log(1 - rn + rns - rs))/((1 - rn + rns - rs)*pow(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs),2))) - 1/((1 - rn + rns - rs)*(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs))),2) + ers2*pow(-(((1/(1 - rs) - 1/(1 - rn + rns - rs))*(log(1 - rn) - log(1 - rn + rns - rs)))/pow(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs),2)) + 1/((1 - rn + rns - rs)*(-log(1 - rn) - log(1 - rs) + log(1 - rn + rns - rs))),2));

        t_k->Fill();
      } // ib
  } // ien

  f_k->cd();
  t_k->Write();
  f_k->Close();
}
