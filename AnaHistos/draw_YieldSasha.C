void draw_YieldSasha(const Int_t process = 0)
{
  gROOT->ProcessLine(".L ReadGraph.C");

  TH1 *h_minv[3][25];  // h_minv[is][ipt]
  for(Int_t is=0; is<3; is++)
    for(Int_t ipt=0; ipt<25; ipt++)
    {
      char name[100];
      sprintf(name,"h_minv_sec%d_pt%d",is,ipt);
      h_minv[is][ipt] = new TH1F(name,name,1000,0.,1.);
    }

  Int_t bin047 = h_minv[0][0]->FindBin(0.047);
  Int_t bin097 = h_minv[0][0]->FindBin(0.097);
  Int_t bin112 = h_minv[0][0]->FindBin(0.112);
  Int_t bin162 = h_minv[0][0]->FindBin(0.162);
  Int_t bin177 = h_minv[0][0]->FindBin(0.177);
  Int_t bin227 = h_minv[0][0]->FindBin(0.227);

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

  TObjArray *Glist = new TObjArray();

  Double_t gx[30], gy[3][30] = {}, egy[3][30] = {};
  Double_t Acceptance[30], eAcceptance[30];
  ReadGraphErrors("Acceptance.root", 3, gx, Acceptance, eAcceptance);

  const Int_t nThread = 10;
  Int_t thread = -1;
  Int_t runnumber;
  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt");

  while( fin >> runnumber )
  {
    thread++;
    if( thread < process*nThread || thread >= (process+1)*nThread ) continue;

    TFile *f = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510MinBias/12233/data/Pi0PP-%d.root",runnumber));
    //TFile *f = new TFile(Form("/phenix/spin/phnxsp01/shura/taxi/Run13pp510MinBias/5096/data/%d.root",runnumber));
    if( f->IsZombie() ) continue;

    TH1 *mchist[3][40];  // mchist[is][ip]
    for(Int_t is=0; is<3; is++)
      for(Int_t ip=0; ip<40; ip++)
      {
        char hname[100];
        sprintf(hname,"mc_s%d_bcc0_pt_%03d_tp",is,5*ip);
        mchist[is][ip] = (TH1*)f->Get(hname);
      }

    for(Int_t is=0; is<3; is++)
      for(Int_t ipt=2; ipt<25; ipt++)
        for(Int_t ip=ptl[ipt]; ip<=ptr[ipt]; ip++)
          h_minv[is][ipt]->Add(mchist[is][ip]);
  }

  for(Int_t is=0; is<3; is++)
    for(Int_t ipt=2; ipt<25; ipt++)
    {
      Double_t npi0 = h_minv[is][ipt]->Integral(bin112,bin162) - ( h_minv[is][ipt]->Integral(bin047,bin097) + h_minv[is][ipt]->Integral(bin177,bin227) ) / 2.;
      npi0 /= 2.;
      gy[is][ipt] = npi0;
      egy[is][ipt] = sqrt(npi0);
    }

  for(Int_t part=0; part<3; part++)
  {
    TGraphErrors *gr = new TGraphErrors(30, gx, gy[part], 0, egy[part]);
    Glist->AddAtAndExpand(gr,part);
    gr->SetName(Form("gr_%d",part));
    gr->SetTitle("#pi^{0} yield");
  }


  mc();
  mcd();
  legi(0, 0.4,0.7,0.7,0.9);

  const char *legname[3] = {"PbScW", "PbScE", "PbGlE"};
  for(Int_t part=0; part<3; part++)
  {
    TGraphErrors *gr = (TGraphErrors*)Glist->At(part);
    aset(gr, "p_{T} [GeV]","Yield", 0.,20., 1.,1e5);
    gPad->SetLogy();
    style(gr, part+20, part+1);
    if(part == 0)
      gr->Draw("AP");
    else
      gr->Draw("P");
    leg0->AddEntry(gr, legname[part], "P");
  }
  leg0->Draw();
  c0->Print("YieldSasha.pdf");

  TFile *fout = new TFile("YieldSasha.root", "RECREATE");
  Glist->Write();
  fout->Close();
}
