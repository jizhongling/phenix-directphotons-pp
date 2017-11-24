void draw_Theta_CV()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos-ertc-cv/total.root");
  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");

  TAxis *axis_sec = (TAxis*)hn_2photon->GetAxis(0);
  TAxis *axis_e = (TAxis*)hn_2photon->GetAxis(2);
  TAxis *axis_th = (TAxis*)hn_2photon->GetAxis(3);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2400);
  gStyle->SetOptStat(0);
  c->Divide(4,4);

  //axis_sec->SetRange(1,6);
  axis_sec->SetRange(7,8);
  
  Int_t ipad = 1;
  Int_t icol = 1;

  for(Int_t iE=0; iE<10; iE++)
  {
    c->cd(ipad++);
    icol = 1;
    bool first = true;

    axis_e->SetRange(iE+1,iE+1);
    Double_t clusterE_low = axis_e->GetBinLowEdge(iE+1);
    Double_t clusterE_high = axis_e->GetBinLowEdge(iE+2);

    for(Int_t ith=0; ith<10; ith++)
    {
      axis_th->SetRange(ith+1,ith+1);
      TH1D *h_2photon = (TH1D*)hn_2photon->Projection(1);
      Int_t nEntries = h_2photon->GetEntries();
      if(nEntries <= 0.)
      {
        h_2photon->Delete();
        icol++;
        continue;
      }
      h_2photon->Scale(1./nEntries);
      cout << "iE=" << iE << "\tith=" << ith << "\tn=" << nEntries << endl;
      h_2photon->SetLineColor(icol);
      if(first)
      {
        first = false;
        char title[100];
        sprintf(title, "PbSc_ClusterE_%4.2f_%4.2f", clusterE_low, clusterE_high);
        h_2photon->SetTitle(title);
        h_2photon->DrawCopy();
      }
      else
      {
        h_2photon->DrawCopy("SAME");
      }
      h_2photon->Delete();
      icol++;
    }
  }

  c->cd(ipad++);
  icol = 1;
  Double_t y = 0.9;
  for(Int_t ith=0; ith<10; ith++)
  {
    Double_t theta_cv = axis_th->GetBinLowEdge(ith+1);
    char buf[100];
    sprintf(buf, "%5.3f", theta_cv);

    TLine *line = new TLine(0.2, y, 0.4, y);
    line->SetLineColor(icol);
    line->Draw();

    TLatex *t = new TLatex();
    t->SetTextFont(22);
    t->SetTextAlign(12);
    t->SetNDC();
    t->DrawLatex(0.5, y, buf);

    y -= 0.05;
    icol++;
  }

  //c->Print("plots/Theta_CV_PbSc.pdf");
  c->Print("plots/Theta_CV_PbGl.pdf");
}
