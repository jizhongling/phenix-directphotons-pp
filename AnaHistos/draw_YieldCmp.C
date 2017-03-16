void draw_YieldCmp()
{
  gROOT->ProcessLine(".L ReadGraph.C");
  const Int_t trig = 2;

  const Double_t gx[30];
  Double_t gy[3][30];
  Double_t egy[3][30];
  Double_t gy2[3][30];
  Double_t egy2[3][30];
  for(Int_t part=0; part<3; part++)
  {
    ReadGraphErrors(Form("Yield-ert%c.root",97+trig), 12*trig+8+part, gx, (Double_t*)gy[part], (Double_t*)egy[part]);
    ReadGraphErrors(Form("Yield-ert%c-sasha.root",97+trig), 12*trig+8+part, gx, (Double_t*)gy2[part], (Double_t*)egy2[part]);
  }

  Double_t Sasha[3][30]= {
    {  0, 0, 105654,  465567, 1890252, 3411618, 3593735, 2809537, 1911005, 1213384,  749961,  464567,  288484,  183085,  119747,   79725,   53802,   38088,   26991,   19584,   39813,   12543,    3879,    1344,     440, 117, 31},
    {   0, 0, 53562,  217808,  843895, 1618568, 1765047, 1431740, 1007771,  660975,  419151,  263268,  166118,  106871,   69846,   46500,   32378,   22462,   16062,   11462,   23313,    7493,    2442,     810,     235,  66, 11},
    {  0, 0,  2389,    7429,   26013,   88324,  162390,  200655,  195866,  165619,  127711,   93753,   67730,   47929,   34037,   24255,   17612,   12492,    9021,    6999,   14566,    5294,    2176,    1096,     509, 225, 100, 62, 16, 10}
  };

  TGraphErrors *gr_ratio[3];
  Double_t rgy[3][30] = {};
  Double_t ergy[3][30] = {};

  for(Int_t part=0; part<3; part++)
  {
    for(Int_t ipt=0; ipt<30; ipt++)
    {
      //rgy[part][ipt] = ( gy[part][ipt] - Sasha[part][ipt] ) / ( Sasha[part][ipt] + 1e-9 );
      //ergy[part][ipt] = egy[part][ipt] / ( Sasha[part][ipt] + 1e-9 );
      rgy[part][ipt] = ( gy[part][ipt] - gy2[part][ipt]/2. ) / ( gy2[part][ipt]/2. + 1e-9 );
      ergy[part][ipt] = gy[part][ipt] / ( gy2[part][ipt]/2. +1e-9 ) * sqrt( pow(egy[part][ipt]/gy[part][ipt],2.) + pow(egy2[part][ipt]/gy2[part][ipt],2.) );
    }
    gr_ratio[part] =  new TGraphErrors(30, gx, rgy[part], 0, ergy[part]);
  }

  TCanvas *c = new TCanvas("c", "Canvas", 600, 600);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
  const char *pname[3] = {"PbScW", "PbScE", "PbGlE"};
  for(Int_t part=0; part<3; part++)
  {
    gr_ratio[part]->SetTitle("Ratio;p_{T} [GeV];ratio;");
    gr_ratio[part]->GetXaxis()->SetRangeUser(1., 20.);
    gr_ratio[part]->GetYaxis()->SetRangeUser(-1., 1.2);
    gr_ratio[part]->SetMarkerStyle(20+part);
    gr_ratio[part]->SetMarkerColor(1+part);
    if(part==0)
      gr_ratio[part]->Draw("APE");
    else
      gr_ratio[part]->Draw("PE");
    leg->AddEntry(gr_ratio[part], Form("%s",pname[part]), "P");
  }
  //leg->Draw();

  c->Print(Form("YieldCmp-ert%c.pdf",97+trig));
}
