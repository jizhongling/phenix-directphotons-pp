#include <AnaToolsTowerID.h>

void draw_Warnmap()
{
  TTree *t1 = new TTree("t1", "Sasha warnmap");

  unsigned int nWarnmap = 0;
  unsigned int nBadSc = 0;
  unsigned int nBadGl = 0;

  int ich = 0;
  int sector = 0;
  int biny = 0;
  int binz = 0;
  int status = 0;

  t1->Branch("sector", &sector, "sector/I");
  t1->Branch("biny", &biny, "biny/I");
  t1->Branch("binz", &binz, "binz/I");
  t1->Branch("status", &status, "status/I");

  ifstream fin("/phenix/plhf/zji/install/share/DirectPhotonPP/warn_all_run13pp500gev.dat");

  while( fin >> ich >> status )
  {
    /* Attention!! I use my indexing for warn map in this program!!! */
    if( ich >= 10368 && ich < 15552 ) { // PbSc
      if( ich < 12960 ) ich += 2592;
      else              ich -= 2592;
    }
    else if( ich >= 15552 )           { // PbGl
      if( ich < 20160 ) ich += 4608;
      else              ich -= 4608;
    }

    /* Get tower location */
    TowerLocation(ich, sector, biny, binz);

    /* Count tower with bad status for PbSc and PbGl */
    if ( status > 0 )
    {
      if( sector < 6 ) nBadSc++;
      else nBadGl++;
    }
    nWarnmap++;

    t1->Fill();
  }
  cout << "NBad PbSc: " << nBadSc << ", PbGl: " << nBadGl
    << ", Total: " << nWarnmap << endl;

  for(int st=0; st<3; st++)
  {
    mc(st, 4,2);
    for(int sec=0; sec<8; sec++)
    {
      mcd(st, sec+1);
      if(st < 2)
        t1->Draw("biny:binz", Form("sector==%d && status==%d",sec,st), "colz");
      else
        t1->Draw("biny:binz", Form("sector==%d && status>=%d",sec,st), "colz");
    }
  }

  TFile *f_out = new TFile("data/Warnmap.root", "RECREATE");
  for(int st=0; st<3; st++)
    mcw( st, Form("status-%d",st) );
  f_out->Close();
}
