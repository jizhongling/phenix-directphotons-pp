void draw_SSkFactor()
{
  const char *name[2] = {"k_{N}", "k_{S}"};

  TFile *f_k = new TFile("data/SSkFactor.root");
  TTree *t_k = (TTree*)f_k->Get("T");
  double rns, kn, ks, ekn, eks;
  t_k->SetBranchAddress("Rate_true", &rns);
  t_k->SetBranchAddress("kN", &kn);
  t_k->SetBranchAddress("kS", &ks);
  t_k->SetBranchAddress("ekN", &ekn);
  t_k->SetBranchAddress("ekS", &eks);

  int nentries = t_k->GetEntries();
  TGraphErrors *gr_k[2];
  int igr_k = 0;
  for(int ik=0; ik<2; ik++)
    gr_k[ik] = new TGraphErrors(nentries);

  for(int ien=0; ien<nentries; ien++)
  {
    t_k->GetEntry(ien);
    if( TMath::Finite(rns+kn+ks+ekn+eks) &&
        rns > 0.1 && ekn < 0.01 && eks < 0.01 )
    {
      gr_k[0]->SetPoint(igr_k, rns, kn);
      gr_k[1]->SetPoint(igr_k, rns, ks);
      gr_k[0]->SetPointError(igr_k, 0., ekn);
      gr_k[1]->SetPointError(igr_k, 0., eks);
      igr_k++;
    }
  }

  mc(0, 2,1);
  for(int ik=0; ik<2; ik++)
  {
    mcd(0, ik+1);
    gr_k[ik]->Set(igr_k);
    gr_k[ik]->SetTitle(name[ik]);
    aset(gr_k[ik], "BBC_{coin}/CLOCK",name[ik], 0.,1., 0.,0.3);
    gr_k[ik]->Draw("AP");
    gr_k[ik]->Fit("pol0", "Q","", 0.,1.);
  }
  c0->Print("plots/SSkFactor.pdf");
}
