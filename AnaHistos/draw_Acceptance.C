TGraphErrors *DivideGraph(TGraphErrors *g1, TGraphErrors *g2)
{
  Int_t gn = g1->GetN();
  Double_t *gx = new Double_t[gn];
  Double_t *gy = new Double_t[gn];
  Double_t *egy = new Double_t[gn];

  for(Int_t i=0; i<gn; i++)
  {
    gx[i] = gy[i] = egy[i] = 0.;

    Double_t g1x, g1y, eg1y;
    Double_t g2x, g2y, eg2y;
    g1->GetPoint(i, g1x, g1y);
    eg1y = g1->GetErrorY(i);
    g2->GetPoint(i, g2x, g2y);
    eg2y = g2->GetErrorY(i);

    if( g1y > 0. && g2y > 0. )
    {
      gx[i] = g1x;
      gy[i] = g1y / g2y;
      egy[i] = gy[i] * sqrt( pow(eg1y/g1y,2.) + pow(eg2y/g2y,2.) );
    }
  }

  TGraphErrors *graph = new TGraphErrors(gn, gx, gy, 0, egy);
  return graph;
}

void GenerateAcceptance(TFile *fsig, TFile *ftot, Int_t ispion)
{
  TH1::SetDefaultSumw2();

  TH1 *h_sig[3];
  TH1 *h_tot;
  if(ispion == 0)
  {
    THnSparse *hn_sig = (THnSparse*)fsig->Get("hn_photon");
    hn_sig->GetAxis(2)->SetRange(1,4);
    h_sig[0] = (TH1*)hn_sig->Projection(0)->Clone("h_sig_PbScW");
    hn_sig->GetAxis(2)->SetRange(5,6);
    h_sig[1] = (TH1*)hn_sig->Projection(0)->Clone("h_sig_PbScE");
    hn_sig->GetAxis(2)->SetRange(7,8);
    h_sig[2] = (TH1*)hn_sig->Projection(0)->Clone("h_sig_PbGlE");
    h_tot = (TH1*)ftot->Get("h_photon");
  }
  else if(ispion == 1)
  {
    THnSparse *hn_sig = (THnSparse*)fsig->Get("hn_photon");
    hn_sig->GetAxis(2)->SetRange(1,4);
    h_sig[0] = (TH1*)hn_sig->Projection(1)->Clone("h_sig_PbScW");
    hn_sig->GetAxis(2)->SetRange(5,6);
    h_sig[1] = (TH1*)hn_sig->Projection(1)->Clone("h_sig_PbScE");
    hn_sig->GetAxis(2)->SetRange(7,8);
    h_sig[2] = (TH1*)hn_sig->Projection(1)->Clone("h_sig_PbGlE");
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
      { gy[i][j] = y_sig / y_tot; egy[i][j] = gy[i][j] * sqrt( pow(ey_sig/y_sig,2.) + pow(ey_tot/y_tot,2.) );
      }
    }
  }

  if(ispion == 0)
  {
    Double_t gy_acc[3][30] = { {0.917729, 1.00633, 1.14825, 1.12254, 1.16406, 1.15534, 1.14387, 1.08724, 1.04397, 1.19259, 1.11401, 1.08243, 1.09937, 1.08396, 1.05962, 1.1155, 1.01567, 1.11697, 1.05258, 1.03892, 1.08092, 1.05828, 1.06569, 1.06397, 1.06614, 1.0498, 1.05094, 1.0626, 1.09115, 1.48079},
      {0.953122, 0.961287, 1.16271, 1.22744, 1.18522, 1.03922, 1.06557, 1.12872, 1.12135, 1.06531, 1.13274, 1.13029, 1.08191, 1.06088, 1.06006, 0.973574, 1.08304, 1.06767, 0.989499, 1.09119, 1.08502, 1.07769, 1.04385, 1.069, 1.05547, 1.03861, 1.09121, 1.05725, 1.04311, 1.5111},
      {1.00682, 0.93535, 1.12194, 1.19194, 1.13168, 1.13583, 1.24479, 1.12686, 1.13593, 1.20425, 1.15742, 1.20454, 1.22982, 1.12088, 1.20276, 1.30014, 1.28389, 1.20988, 1.16998, 1.22525, 1.26973, 1.29905, 1.28267, 1.32714, 1.32478, 1.35913, 1.42392, 1.45609, 1.84625, 3.00666} };
    Double_t egy_acc[3][30] = { {0.0353803, 0.037015, 0.0498044, 0.0507489, 0.0532061, 0.0507488, 0.0510535, 0.0479079, 0.0473473, 0.0523246, 0.0511518, 0.0486528, 0.0497655, 0.0488145, 0.0474966, 0.0501719, 0.0456049, 0.050353, 0.0470434, 0.0472703, 0.0255918, 0.0251296, 0.0259282, 0.0259721, 0.0269765, 0.0271598, 0.0293824, 0.0334561, 0.0431621, 0.0952064},
      {0.0542788, 0.0523834, 0.0756496, 0.080098, 0.0749419, 0.0676348, 0.0711552, 0.0723044, 0.0716856, 0.0672427, 0.0741465, 0.0748741, 0.0704224, 0.0712935, 0.0717517, 0.0664484, 0.0722228, 0.0690087, 0.0654933, 0.0732659, 0.0379495, 0.0374642, 0.0362525, 0.0381689, 0.0384977, 0.0393928, 0.0447302, 0.0488738, 0.0596545, 0.142828},
      {0.0535118, 0.0465713, 0.0677599, 0.0743985, 0.0689997, 0.0709085, 0.0744999, 0.0684519, 0.068995, 0.0726542, 0.0717698, 0.0753872, 0.0765055, 0.0686541, 0.0737477, 0.0801868, 0.0802004, 0.0746237, 0.0710075, 0.0762464, 0.0413532, 0.0426517, 0.0423268, 0.044556, 0.0449116, 0.0483652, 0.0532628, 0.0610961, 0.092239, 0.224152} };
  }
  else if(ispion == 1)
  {
    Double_t gy_acc[3][30] = { {0, 0, 0.88385, 1.04499, 1.0231, 1.13176, 1.07586, 1.04761, 1.087, 1.05929, 1.08917, 1.11467, 1.00637, 1.1191, 0.994612, 1.09007, 1.08234, 1.00524, 1.08042, 1.05133, 1.04236, 1.0592, 1.07974, 1.04979, 1.05967, 1.0518, 1.02682, 1.0698, 1.03018, 0.91226},
      {0, 0, 0.855307, 0.993325, 0.909117, 1.07722, 1.0727, 1.09681, 0.959494, 1.23678, 1.0138, 1.00096, 1.09915, 1.04102, 1.05418, 1.08455, 1.04873, 1.08632, 1.04557, 1.06438, 1.05826, 1.06538, 1.0469, 1.07063, 1.05238, 1.07157, 1.04929, 1.07831, 1.02589, 0.888965},
      {0, 0, 0.826921, 0.993274, 0.961897, 1.01002, 1.06273, 0.97499, 1.04529, 1.11787, 1.04805, 1.10157, 1.0097, 1.21474, 1.14056, 0.978027, 1.16989, 1.17928, 1.17357, 1.14793, 1.11843, 1.19968, 1.20214, 1.21795, 1.22465, 1.27215, 1.26836, 1.28446, 1.23853, 1.27468} };
    Double_t egy_acc[3][30] = { {0, 0, 0.056554, 0.06223, 0.0546813, 0.0574719, 0.0527679, 0.0512593, 0.0508041, 0.0495695, 0.051302, 0.0517243, 0.0467729, 0.0503892, 0.0442187, 0.0483369, 0.0471697, 0.0429209, 0.046725, 0.0451935, 0.0230536, 0.0225947, 0.0227181, 0.0215629, 0.0214289, 0.0212466, 0.0202437, 0.0212453, 0.0199182, 0.0188465},
      {0, 0, 0.0892606, 0.0883507, 0.0723156, 0.0836635, 0.0851724, 0.0802716, 0.0700661, 0.0878699, 0.0671981, 0.0664275, 0.0722757, 0.0675377, 0.0673299, 0.0685775, 0.0668313, 0.0700965, 0.0662028, 0.0660608, 0.0338414, 0.0326503, 0.0320392, 0.0319711, 0.0305577, 0.0311139, 0.030268, 0.0307963, 0.0290396, 0.0262014},
      {0, 0, 0.0753775, 0.0791471, 0.0690368, 0.0726844, 0.0699909, 0.0615804, 0.0642381, 0.0683563, 0.0655393, 0.0671857, 0.0623245, 0.0713398, 0.0665938, 0.0579418, 0.0675888, 0.0672193, 0.0690177, 0.0652343, 0.0330419, 0.0348125, 0.0336409, 0.0341586, 0.0343425, 0.0346996, 0.0345377, 0.0349239, 0.0332483, 0.033976} };
  }

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);

  TGraphErrors *gr1[3];
  TGraphErrors *gr2[3];
  TGraphErrors *gr[3];
  for(Int_t i=0; i<3; i++)
  {
    gr[i] = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
    //gr2[i] = new TGraphErrors(gn, gx[i], gy_acc[i], 0, egy_acc[i]);
    //gr[i] = DivideGraph(gr1[i], gr2[i]);
    if(ispion == 0)
      gr[i]->SetTitle("Photon acceptance");
    else if(ispion == 1)
      gr[i]->SetTitle("#pi^{0} acceptance");
    gr[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr[i]->GetYaxis()->SetTitle("acceptance");
    gr[i]->GetYaxis()->SetTitleOffset(1.2);
    gr[i]->GetXaxis()->SetRangeUser(0., 30.);
    gr[i]->GetYaxis()->SetRangeUser(0., 0.15);
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
    cout << "\nPart " << i << endl;
    cout << "      const Double_t Acceptance[30] = {";
    for(Int_t j=0; j<gn; j++)
      cout << gy[i][j] << ",";
    cout << "};\n";
    cout << "      const Double_t eAcceptance[30] = {";
    for(Int_t j=0; j<gn; j++)
      cout << egy[i][j] << ",";
    cout << "};\n";
  }

  return;
}

void draw_Acceptance()
{
  TFile *fsig = new TFile("MissingRatio-histo.root");
  TFile *ftot = new TFile("AnaPHPythia-histo.root");
  for(Int_t ispion=0; ispion<2; ispion++)
  {
    cout << "\nispion " << ispion << endl;
    GenerateAcceptance(fsig, ftot, ispion);
  }
}
