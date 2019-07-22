void draw_DCDeadmap(int runnumber)
{
  //mc();
  //mcd();
  //gGeoManager = new TGeoManager();
  //TGeoPolygon *pl3 = new TGeoPolygon(3);
  //double x3[3] = {20., 40., 60.};
  //double y3[3] = {-0.5, 0.5, -0.5};
  //pl3->SetXY(x3, y3);
  //pl3->FinishPolygon();
  //pl3->Draw();
  //double xy1[2] = {40., 0.};
  //double xy2[2] = {50., 0.5};
  //cout << "Contains or not: " << pl3->Contains(xy1) << ", " << pl3->Contains(xy2) << endl;
  //return;

  TFile *f_dc = new TFile("data/DCCheck.root");
  TLine *line = new TLine();
  line->SetLineColor(kRed);
  line->SetLineWidth(5);
  line->SetLineStyle(2);

  for(int ns=0; ns<2; ns++)
    for(int we=0; we<2; we++)
    {
      mc(ns*2+we);
      mcd(ns*2+we);
      TString NS = ns ? "S" : "N";
      TString WE = we ? "E" : "W";
      TH2 *h2_board_run = (TH2*)f_dc->Get( Form("h2_board_ns%d_we%d_%d",ns,we,runnumber) );
      if(!h2_board_run) continue;
      h2_board_run->SetTitle( Form("%d, %s",runnumber,(NS+WE).Data()) );
      h2_board_run->DrawCopy("COLZ");
      line->DrawLine(40., -0.5, 40., 0.5);
      delete h2_board_run;
    }
}

// par[4]: x1, y1, x2, y2
double myline(double *x, double *par)
{
  if(TMath::Abs(par[0]-par[2]) < 1e-9)
  {
    double f = (par[1]+par[3])/2.;
    return f;
  }
  double f = (par[3]-par[1])/(par[2]-par[0])*(x[0]-par[0]) + par[1];
  return f;
}

void cutLine(TLine *line, double alphacut1, double alphacut2)
{
  double board1 = line->GetX1();
  double alpha1 = line->GetY1();
  double board2 = line->GetX2();
  double alpha2 = line->GetY2();
  TF1 *f_line = new TF1("f_line", myline, -10., 10., 4);
  f_line->SetParameters(alpha1, board1, alpha2, board2);
  double boardcut1 = f_line->Eval(alphacut1);
  double boardcut2 = f_line->Eval(alphacut2);
  cout.precision(4);
  cout << "(" << boardcut1 << "," << alphacut1 << ")->(" << boardcut2 << "," << alphacut2 << ")" << endl;
}

void getLine(TLine *line)
{
  double board1 = line->GetX1();
  double alpha1 = line->GetY1();
  double board2 = line->GetX2();
  double alpha2 = line->GetY2();
  if(TMath::Abs(alpha1-alpha2) < 1e-9)
  {
    cout << "Too narrow alpha range!" << endl;
    return;
  }
  double k = (board2-board1)/(alpha2-alpha1);
  double b = board1 - k*alpha1;
  cout.precision(4);
  cout << "board = " << k << "*alpha + " << b << endl;
}
