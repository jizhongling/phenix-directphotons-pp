// mysrc64 new
// g++ -std=c++11 -Wall -I$MYINSTALL/include -L$MYINSTALL/lib `root-config --cflags --glibs` -lLHAPDF -o anaPDF_xdg anaPDF_xdg.cc
#include <cmath>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include "LHAPDF/LHAPDF.h"

using namespace std;

int main()
{
  const char *pdfname[3] = {"DSSV", "JAM with W", "JAM without W"};
  const char *pdfset[4] = {"DSSV_REP_LHAPDF6", "JAM22ppdf", "JAM22_pol_SU23_pos_g", "JAM22_pol_SU23_neg_g"};

  auto c0 = new TCanvas("c0", "c0", 3*600, 600);
  c0->Divide(3, 1);

  for(int iset=0; iset<4; iset++)
  {
    c0->cd(iset<3?iset+1:3);
    vector<LHAPDF::PDF*> v_pdf = LHAPDF::mkPDFs(pdfset[iset]);
    const int nrep = stoi(v_pdf.at(0)->info().get_entry("NumMembers")) - 1;
    for(int irep=1; irep<=nrep; irep++)
    {
      auto gr_xg = new TGraph(101);
      for(int ix=0; ix<101; ix++)
      {
        double log10x = -3. + 0.03 * (double)ix;
        double x = pow(10., log10x);
        double xg = v_pdf.at(irep)->xfxQ2(21, x, 10.);
        gr_xg->SetPoint(ix, log10x, xg);
      } // ix
      gr_xg->SetTitle(pdfname[iset]);
      gr_xg->GetXaxis()->SetTitle("log10(x)");
      gr_xg->GetYaxis()->SetTitle("x#Deltag");
      gr_xg->GetXaxis()->SetRangeUser(-3., 0.);
      gr_xg->GetYaxis()->SetRangeUser(-0.5, 0.5);
      gr_xg->Draw(iset<3&&irep==1?"AC":"C");
    } // irep
  } // iset

  c0->Print("plots/PDF-xdg.pdf");

  return 0;
}
