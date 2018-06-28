void draw_Acceptance_Cmp()
{
  const char *fname[4] = {"PH-nowarn", "PH-warn", "Fast-nowarn", "Fast-warn"};
  const char *hname[4] = {"input", "geom", "goodtower", "cut"};

  mc(0, 2,4);
  legi(0, 0.2,0.2,0.9,0.9);

  for(int i=0; i<4; i++)
  {
    TFile *f = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/AnaFastMC-%s-histo.root",fname[i]));

    THnSparse *hn_pion[4];

    for(int ih=0; ih<4; ih++)
    {
      hn_pion[ih] = (THnSparse*)f->Get(Form("hn_pion_%s",hname[ih]));

      if( i==0 && ih<3 )
      {
        for(int ip=0; ip<2; ip++)
        {
          TH1 *h_proj = hn_pion[ih]->Projection(ip+1);
          aset(h_proj);
          h_proj->SetMinimum(0.);
          h_proj->SetLineColor(3*i+ih+1);
          if(ih==0)
          {
            mcd(0, ip+1);
            h_proj->SetTitle("Input, geom, goodtower without warnmap");
            h_proj->Draw("HIST");
          }
          if(ih!=0)
          {
            mcd(0, ip+1);
            h_proj->Draw("HIST SAME");
          }
          if(ih==2)
          {
            mcd(0, ip+3);
            h_proj->SetTitle("Goodtower without warnmap, and goodtower and additional cuts with warnmap");
            h_proj->Draw("HIST");
          }
        }
        leg0->AddEntry(h_proj, Form("%s-%s",fname[i],hname[ih]));
      }

      if( i==1 && ih>1 )
      {
        for(int ip=0; ip<2; ip++)
        {
          mcd(0, ip+3);
          TH1 *h_proj = hn_pion[ih]->Projection(ip+1);
          aset(h_proj);
          h_proj->SetMinimum(0.);
          h_proj->SetLineColor(3*i+ih+1);
          h_proj->Draw("HIST SAME");
        }
        leg0->AddEntry(h_proj, Form("%s-%s",fname[i],hname[ih]));
      }

      if( ( i==0 || i==3 ) && ih == 0 )
      {
        for(int ip=0; ip<2; ip++)
        {
          mcd(0, ip+5);
          TH1 *h_proj = hn_pion[ih]->Projection(ip+1);
          aset(h_proj);
          h_proj->SetMinimum(0.);
          h_proj->SetMaximum(18000.);
          h_proj->SetLineColor(i+ih+1);
          h_proj->SetTitle("Total input for PHParticleGen and FastMC");
          if(i==0)
            h_proj->Draw("HIST");
          else
            h_proj->Draw("HIST SAME");
        }
        if(i==3)
          leg0->AddEntry(h_proj, Form("%s-%s",fname[i],hname[ih]));
      }
    }
  }
  
  mcd(0, 7);
  leg0->Draw();

  c0->Print("plots/Acceptance-cmp.pdf");
}
