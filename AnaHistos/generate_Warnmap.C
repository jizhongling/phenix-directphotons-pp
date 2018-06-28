void generate_Warnmap()
{
  TFile *f = new TFile("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/AnalysisTrain/Run13_EMC_TOF_Calibration_YIS/Run13pp510_WarnMap_05.root");
  TTree *T = (TTree*)f->Get("T");

  int warnmap[8][48][96];
  T->SetBranchAddress("warnmap", warnmap);
  T->GetEntry(0);

  ofstream fout("/phenix/plhf/zji/install/share/DirectPhotonPP/Warnmap_Run13pp510_ISYOON.txt", ofstream::trunc);

  for(int i=0; i<8; i++)
  {
    int sector = i<4 ? i : 11-i;
    for(int iy=0; iy<48; iy++)
      for(int iz=0; iz<96; iz++)
      {
        if(sector<6 && (iy>=36 || iz>=72) )
          continue;

        int status = 150;
        // Hot
        if(warnmap[i][iy][iz] & 1)
          status = 50;
        // Dead
        else if(warnmap[i][iy][iz] & 2)
          status = 100;
        // Edge
        else if(warnmap[i][iy][iz] & 4)
          status = 10;
        // Neighbor
        else if(warnmap[i][iy][iz] & 8)
          status = 40;
        // Good
        else if(!warnmap[i][iy][iz])
          status = 0;

        fout << sector << " " << iy << " " << iz << " " << status << endl;
      }
  }

  fout.close();
}
