// mysrc64 new
// g++ -std=c++11 -Wall -I$MYINSTALL/include -L$MYINSTALL/lib `root-config --cflags --glibs` -lLHAPDF -o anaPDF_xdg anaPDF_xdg.cc
#include <cmath>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TAxis.h>
#include "LHAPDF/LHAPDF.h"

using namespace std;

int main()
{
  const int npdf = 4;
  const int nset = 5;
  const char *pdfname[npdf] = {"DSSV14", "JAM22 with W", "JAM22 positivity without W", "JAM22 without W"};
  const char *pdfset[nset] = {"DSSV_REP_LHAPDF6", "JAM22ppdf", "JAM22_pol_positivity", "JAM22_pol_SU23_pos_g", "JAM22_pol_SU23_neg_g"};

  auto c0 = new TCanvas("c0", "c0", 2*600, 2*600);
  c0->Divide(2, 2);

  for(int iset=0; iset<nset; iset++)
  {
    c0->cd(iset<npdf?iset+1:npdf);
    gPad->SetLogx();
    vector<LHAPDF::PDF*> v_pdf = LHAPDF::mkPDFs(pdfset[iset]);
    const int irep_start = iset == 0 ? 1 : 0;
    const int irep_end = stoi(v_pdf.at(0)->info().get_entry("NumMembers")) - 1;
    for(int irep=irep_start; irep<=irep_end; irep++)
    {
      auto gr_xg = new TGraph(101);
      for(int ix=0; ix<101; ix++)
      {
        double log10x = -3. + 0.03 * (double)ix;
        double x = pow(10., log10x);
        double xg = v_pdf.at(irep)->xfxQ2(21, x, 10.);
        gr_xg->SetPoint(ix, x, xg);
      } // ix
      gr_xg->SetTitle(pdfname[iset]);
      gr_xg->GetXaxis()->SetTitle("x");
      gr_xg->GetYaxis()->SetTitle("x#Deltag");
      gr_xg->GetXaxis()->SetRangeUser(1e-3, 1.);
      gr_xg->GetYaxis()->SetRangeUser(-0.5, 0.5);
      gr_xg->Draw(iset<npdf&&irep==irep_start?"AC":"C");
    } // irep
  } // iset

  c0->Print("plots/PDF-xdg.pdf");

  return 0;
}
