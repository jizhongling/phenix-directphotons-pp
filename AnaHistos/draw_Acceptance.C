void GenerateAcceptance(TFile *fsig, TFile *ftot, Int_t ispion)
{
  TH1::SetDefaultSumw2();

  TH1 *h_sig[3];
  TH1 *h_tot;
  if(ispion == 0)
  {
    THnSparse *hn_sig = (THnSparse*)fsig->Get("hn_photon");
    hn_sig->GetAxis(1)->SetRange(1,4);
    h_sig[0] = (TH1*)hn_sig->Projection(0)->Clone("h_sig_PbScW");
    hn_sig->GetAxis(1)->SetRange(5,6);
    h_sig[1] = (TH1*)hn_sig->Projection(0)->Clone("h_sig_PbScE");
    hn_sig->GetAxis(1)->SetRange(7,8);
    h_sig[2] = (TH1*)hn_sig->Projection(0)->Clone("h_sig_PbGlE");
    h_tot = (TH1*)ftot->Get("h_photon");
  }
  else if(ispion == 1)
  {
    TH2 *h2_sig = (TH2*)fsig->Get("h2_pion");
    h_sig[0] = (TH1*)h2_sig->ProjectionX("_px",1,4)->Clone("h_sig_PbScW");
    h_sig[1] = (TH1*)h2_sig->ProjectionX("_px",5,6)->Clone("h_sig_PbScE");
    h_sig[2] = (TH1*)h2_sig->ProjectionX("_px",7,8)->Clone("h_sig_PbGlE");
    h_tot = (TH1*)ftot->Get("h_pion");
  }

  Int_t gn = 30;
  Double_t *gx[3];
  Double_t *gy[3];
  Double_t *egy[3];
  for(Int_t i=0; i<3; i++)
  {
    gx[i] = new Double_t[gn];
    gy[i] = new Double_t[gn];
    egy[i] = new Double_t[gn];
    for(Int_t j=0; j<30; j++)
    {
      gx[i][j] = gy[i][j] = egy[i][j] = 0.;
      gx[i][j] = h_tot->GetXaxis()->GetBinCenter(j+1);
      y_sig = h_sig[i]->GetBinContent(j+1);
      y_tot = h_tot->GetBinContent(j+1);
      ey_sig = h_sig[i]->GetBinError(j+1);
      ey_tot = h_tot->GetBinError(j+1);
      if( y_sig > 0. && y_tot > 0.)
      {
        gy[i][j] = y_sig / y_tot;
        egy[i][j] = gy[i][j] * sqrt( pow(ey_sig/y_sig,2.) + pow(ey_tot/y_tot,2.) );
      }
    }
  }

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TGraphErrors *gr[3];
  for(Int_t i=0; i<3; i++)
  {
    gr[i] = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
    if(ispion == 0)
      gr[i]->SetTitle("Photon acceptance");
    else if(ispion == 1)
      gr[i]->SetTitle("#pi^{0} acceptance");
    gr[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr[i]->GetYaxis()->SetTitle("acceptance");
    gr[i]->GetYaxis()->SetTitleOffset(1.2);
    gr[i]->GetXaxis()->SetRangeUser(0., 30.);
    gr[i]->GetYaxis()->SetRangeUser(0., 0.1);
    gr[i]->SetMarkerColor(i+1);
    gr[i]->SetMarkerStyle(i+20);
    if(i==0)
      gr[i]->Draw("AP");
    else
      gr[i]->Draw("P");
  }

  TLegend *leg = new TLegend(0.4, 0.7, 0.6, 0.9);
  leg->AddEntry(gr[0], "PbScW", "LEP");
  leg->AddEntry(gr[1], "PbScE", "LEP");
  leg->AddEntry(gr[2], "PbGlE", "LEP");
  leg->Draw();

  if(ispion == 0)
  {
    c->Print("Acceptance-photon.pdf");
    delete c;
  }
  else if(ispion == 1)
  {
    c->Print("Acceptance-pion.pdf");
    delete c;
  }

  for(Int_t i=0; i<3; i++)
  {
    cout << "\nPart " << i << " Y: ";
    for(Int_t j=0; j<gn; j++)
      cout << gy[i][j] << ",";
  }
  cout << endl;

  for(Int_t i=0; i<3; i++)
  {
    cout << "\nPart " << i << " EY: ";
    for(Int_t j=0; j<gn; j++)
      cout << egy[i][j] << ",";
  }
  cout << endl;

  return;
}

void draw_Acceptance()
{
  TFile *fsig = new TFile("Acceptance-signal.root");
  TFile *ftot = new TFile("Acceptance-total.root");
  for(Int_t ispion=0; ispion<2; ispion++)
  {
    cout << "\nispion " << ispion << endl;
    GenerateAcceptance(fsig, ftot, ispion);
  }
}
