void draw_Pileup()
{
  const char *name[4] = {"PbSc without ToF", "PbGl without ToF", "PbSc with ToF", "PbGl with ToF"};
  TGraphErrors *gr[4];
  TMultiGraph *mg[4];
  for(Int_t img=0; img<4; img++)
    mg[img] = new TMultiGraph();

  for(Int_t i=0; i<45; i++)
  {
    TFile *f = new TFile(Form("pileup/Sasha-%d.root",i));
    if(f->IsZombie()) continue;

    for(Int_t img=0; img<4; img++)
    {
      gr[img] = (TGraphErrors*)f->Get(Form("gr_%d",img));
      mg[img]->Add(gr[img]);
    }
  }

  mc(0, 2,2);
  for(Int_t img=0; img<4; img++)
  {
    mcd(0, img+1);
    mg[img]->Draw("AP");  // must before GetXaxis()
    mg[img]->SetTitle(name[img]);
    mg[img]->GetXaxis()->SetTitle("Nmb/Nclock");
    mg[img]->GetYaxis()->SetTitle("Npi0/Nmb");
    mg[img]->GetXaxis()->SetLimits(0., 0.2);  // Do not use SetRangeUser()
  }
  //mg[0]->GetYaxis()->SetRangeUser(0., 4.5e-3);  // Do not use SetLimits()
  //mg[1]->GetYaxis()->SetRangeUser(0., 0.5e-3);  // Do not use SetLimits()
  mg[2]->GetYaxis()->SetRangeUser(0., 0.8e-3);  // Do not use SetLimits()
  mg[3]->GetYaxis()->SetRangeUser(0., 0.8e-3);  // Do not use SetLimits()
  c0->Print("Pileup-Sasha.pdf");
}
