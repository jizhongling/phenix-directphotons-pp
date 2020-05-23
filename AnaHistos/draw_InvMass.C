void draw_InvMass()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/PhotonHistos-Sasha.root");
  THnSparse *hn_inv_mass_2photon = (THnSparse*)f->Get("inv_mass_2photon");

  TAxis *axis0 = (TAxis*)hn_inv_mass_2photon->GetAxis(0);
  //TAxis *axis1 = (TAxis*)hn_inv_mass_2photon->GetAxis(1);
  TAxis *axis2 = (TAxis*)hn_inv_mass_2photon->GetAxis(2);
  TAxis *axis3 = (TAxis*)hn_inv_mass_2photon->GetAxis(3);

  int Last0 = axis0->GetLast();
  //int Last1 = axis1->GetLast();
  int Last2 = axis2->GetLast();
  int Last3 = axis3->GetLast();

  int nsec = (Last2<8) ? Last2 : 8;
  for(int isec=0; isec<nsec; isec++)
  {
    TCanvas *c = new TCanvas("c", "Canvas", 1200, 1200);
    gStyle->SetOptStat(0);
    c->Divide(4,4);

    axis2->SetRange(isec+1,isec+1);
    int npt = (Last0<15) ? Last0 : 15;
    for(int ipt=1; ipt<=npt; ipt++)
    {
      c->cd(ipt);
      axis0->SetRange(ipt+1,ipt+1);
      double pt_low = axis0->GetBinLowEdge(ipt+1);
      double pt_high = axis0->GetBinLowEdge(ipt+2);
      for(int icon=2; icon<=Last3; icon++)
      {
        axis3->SetRange(icon,icon);
        TH1D *hnp_inv_mass_2photon = (TH1D*)hn_inv_mass_2photon->Projection(1);
        hnp_inv_mass_2photon->SetLineColor(icon);
        if(icon==2)
        {
          char title[100];
          sprintf(title, "Sector_%d_pT_%4.2f_%4.2f", isec, pt_low, pt_high);
          hnp_inv_mass_2photon->SetTitle(title);
          hnp_inv_mass_2photon->DrawCopy();
        }
        else
        {
          hnp_inv_mass_2photon->DrawCopy("SAME");
        }
        hnp_inv_mass_2photon->Delete();
      }
    }

    c->cd(npt+1);
    const char *cond[] = {"Direct Photon", "Photon", "E_{min}", "ToF", "Shape", "#theta_{CV}"};

    double y = 0.8;
    for(int icon=2; icon<=Last3; icon++)
    {

      TLine *line = new TLine(0.2, y, 0.5, y);
      line->SetLineColor(icon);
      line->Draw();

      TLatex *t = new TLatex();
      t->SetTextFont(22);
      t->SetTextAlign(12);
      t->SetNDC();
      t->DrawLatex(0.6, y, cond[icon-1]);

      y -= 0.1;
    }

    char buf[100];
    sprintf(buf, "InvMass-%d.pdf", isec);
    c->Print(buf);
    delete c;
  }
}
