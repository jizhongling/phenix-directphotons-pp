void draw_ALL()
{
  const Int_t npT = 30;
  const Float_t pTbins[npT+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0};

  //Float_t ALLbar[npT] = {};
  //Float_t eALLbar[npT] = {};

  for(Int_t itype=0; itype<3; itype++)
  {
    TCanvas *c = new TCanvas("c", "Canvas", 2000, 2400);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();
    c->Divide(5,6);

    for(Int_t ipt=0; ipt<npT; ipt++)
    {
      c->cd(ipt+1);
      bool drawAxis = true;
      for(Int_t ipart=0; ipart<3; ipart++)
        for(Int_t icr=0; icr<2; icr++)
        {
          TGraphErrors *gr = new TGraphErrors(Form("ALL/type%d-part%d-crossing%d-pT%d.txt",itype,ipart,icr,ipt), "%lg %lg %lg");
          if( gr->GetN() == 0 ) continue;
          gr->SetTitle(Form("p_{T} %.2f-%.2f", pTbins[ipt], pTbins[ipt+1]));
          gr->GetXaxis()->SetTitle("p_{T} [GeV]");
          gr->GetYaxis()->SetTitle("A_{LL}");
          gr->GetYaxis()->SetTitleOffset(1.2);
          gr->GetXaxis()->SetRangeUser(386500, 398500);
          gr->GetYaxis()->SetRangeUser(-5., 5.);
          gr->SetMarkerStyle(20+ipart);
          gr->SetMarkerColor(1+icr);
          //gr->Fit("pol0", "Q");
          //TF1 *fit = gr->GetFunction("pol0");
          //ALLbar[ipt] = fit->GetParameter(0);
          //eALLbar[ipt] = fit->GetParError(0);
          if(drawAxis)
          {
            gr->Draw("AP");
            drawAxis = false;
          }
          else
          {
            gr->Draw("P");
          }
        }
    }

    c->Print(Form("ALL-runbyrun-type%d.pdf",itype));
    delete c;
  }

  //TCanvas *c0 = new TCanvas("c0", "Canvas", 600, 600);
  //gStyle->SetOptStat(0);

  //Float_t xpt[npT];
  //for(Int_t ipt=0; ipt<npT; ipt++)
  //  xpt[ipt] = ( pTbins[ipt] + pTbins[ipt+1] ) / 2.;

  //TGraphErrors *gr0 = new TGraphErrors(npT, xpt, ALLbar, 0, eALLbar);
  //gr0->SetTitle("Direct Photon A_{LL}");
  //gr0->GetXaxis()->SetTitle("p_{T} [GeV]");
  //gr0->GetYaxis()->SetTitle("A_{LL}");
  //gr0->GetYaxis()->SetTitleOffset(1.2);
  //gr0->GetXaxis()->SetRangeUser(0., 15.);
  //gr0->GetYaxis()->SetRangeUser(-0.05, 0.08);
  //gr0->SetMarkerStyle(20);
  //gr0->SetMarkerColor(2);
  //gr0->Draw("AP");

  //c0->Print("ALL-average.pdf");
}
