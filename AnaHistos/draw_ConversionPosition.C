void draw_ConversionPosition()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/MissingRatio-histo.root");
  TH2 *h2_radius= (TH2*)f->Get("h2_radius");
  TH2 *h2_angle= (TH2*)f->Get("h2_angle");
  THnSparse *hn_position = (THnSparse*)f->Get("hn_conversion_position");

  for(int ic=0; ic<4; ic++)
    mc(ic, 6,5);

  int ipad = 1;
  for(int ipt=0; ipt<30; ipt++)
  {
    double pTlow = h2_radius->GetXaxis()->GetBinLowEdge(ipt+1);
    double pTup = h2_radius->GetXaxis()->GetBinUpEdge(ipt+1);

    mcd(0, ipad);
    TH1 *h_radius= h2_radius->ProjectionY(Form("h_radius_py_%d",ipt), ipt+1,ipt+1);
    h_radius->SetTitle(Form("pT: %3.1f-%3.1f GeV",pTlow,pTup));
    aset(h_radius);
    h_radius->DrawCopy("HIST");

    mcd(1, ipad);
    aset(h_radius, "","", -50.,50.);
    h_radius->DrawCopy("HIST");
    TLine *line = new TLine();
    line->SetLineColor(kRed);
    double ymax = h_radius->GetYaxis()->GetXmax();
    ymax = 1000.;
    line->DrawLine(2.63,0., 2.63,ymax);
    line->DrawLine(5.13,0., 5.13,ymax);
    line->DrawLine(11.77,0., 11.77,ymax);
    line->DrawLine(16.69,0., 16.69,ymax);

    mcd(2, ipad);
    TH1 *h_angle = h2_angle->ProjectionY(Form("h_angle_py_%d",ipt), ipt+1,ipt+1);
    h_angle->SetTitle(Form("pT: %3.1f-%3.1f GeV",pTlow,pTup));
    aset(h_angle);
    h_angle->DrawCopy("HIST");

    mcd(3, ipad);
    hn_position->GetAxis(0)->SetRange(ipt+1,ipt+1);
    TH2 *h2_position = hn_position->Projection(2,1);
    h2_position->SetTitle(Form("pT: %3.1f-%3.1f GeV",pTlow,pTup));
    aset(h2_position);
    h2_position->DrawCopy();
    delete h2_position;

    ipad++;
  }

  c0->Print("plots/Conversion-radius-full.pdf");
  c1->Print("plots/Conversion-radius-center.pdf");
  c2->Print("plots/Conversion-angle.pdf");
  c3->Print("plots/Conversion-position.pdf");
}
