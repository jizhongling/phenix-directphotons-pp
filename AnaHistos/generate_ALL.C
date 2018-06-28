void GetRawAsymm(TFile *f, int part, int pattern, Long64_t *yield)
{
  const double BR[30] = {0.0444823,0.155126,0.189877,0.209775,0.227416,0.235346,0.242619,0.258304,0.254693,0.219828,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313};

  if(part == 0)
  {
    const int sector_low = 1;
    const int sector_high = 4;
    const double MissR[30] = {0.576394, 1.5215, 2.09614, 1.55921, 1.98752, 1.47444, 1.54181, 1.73786, 1.36634, 1.40074, 1.27356, 1.14209, 1.11697, 1.03645, 0.886776, 1.00954, 0.928018, 0.729948, 0.753175, 0.634566, 0.928788, 1.14879, 1.18122, 0.905212, 1.30759, 1.91196, 2.2468, 3.98329, 5.55577, 7.05814};
    const double TrigE[30] = {0.000164287,0.000226765,0.000554624,0.00132205,0.00273391,0.00653341,0.0205823,0.0836921,0.247432,0.506061,0.711566,0.897269,0.920601,0.926923,0.952096,0.963303,0.897059,0.916667,0.931034,1,0.930233,0.916667,1,1,1,1,1,1,1,1};
  }
  else if(part == 1)
  {
    const int sector_low = 5;
    const int sector_high = 6;
    const double MissR[30] = {2.58414, 1.4808, 1.60466, 1.45547, 1.53042, 1.18097, 1.23304, 1.54858, 0.916936, 1.45538, 1.00605, 1.38998, 0.851646, 0.760959, 0.78359, 0.605272, 0.398232, 0.43892, 0.425398, 0.421566, 0.740498, 0.681352, 0.782055, 0.615879, 1.18684, 1.0558, 2.03564, 3.30069, 4.03718, 7.15849};
    const double TrigE[30] = {0.000170962,0.000212571,0.000559405,0.00126382,0.00254494,0.00865695,0.0343194,0.112287,0.26078,0.447848,0.645638,0.836788,0.864,0.934641,0.97561,0.941176,0.913043,0.96875,1,0.909091,0.90625,0.888889,1,1,1,1,1,1,1,1};
  }
  else if(part == 2)
  {
    const int sector_low = 7;
    const int sector_high = 8;
    const double MissR[30] = {0.579534, 1.17435, 2.05909, 1.63442, 0.980906, 1.18976, 1.32591, 1.34993, 1.49276, 1.03652, 0.962146, 0.69767, 0.897513, 0.685265, 0.684048, 0.536245, 0.442107, 0.455011, 0.48081, 0.561676, 0.650728, 0.694951, 0.459946, 0.629772, 0.565775, 0.70639, 0.66377, 1.10261, 1.46038, 4.89179};
    const double TrigE[30] = {0.000179314,0.000245431,0.000527023,0.00106538,0.00213491,0.00297094,0.0233425,0.0625663,0.143363,0.257883,0.379648,0.516522,0.674352,0.665179,0.690141,0.75,0.788732,0.790698,0.83871,0.818182,0.829787,0.933333,1,1,1,1,1,1,1,1};
  }
  else
  {
    cerr << "Wrong part " << part << endl;
    exit(1);
  }

  if(pattern == 1)
  {
    const int pattern_bin = 3;
  }
  else if(pattern == -1)
  {
    const int pattern_bin = 1;
  }
  else
  {
    cout << "Irrelevant pattern " << pattern << endl;
    return;
  }

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  hn_1photon->GetAxis(0)->SetRange(sector_low, sector_high);
  hn_1photon->GetAxis(2)->SetRange(pattern_bin, pattern_bin);
  TH1 *h_1photon = hn_1photon->Projection(1);

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");

  TAxis *axis0 = hn_2photon->GetAxis(0);
  TAxis *axis1 = hn_2photon->GetAxis(1);
  TAxis *axis2 = hn_2photon->GetAxis(2);
  TAxis *axis3 = hn_2photon->GetAxis(3);

  axis0->SetRange(sector_low, sector_high);
  axis3->SetRange(pattern_bin, pattern_bin);

  int bin047 = axis2->FindBin(0.047);
  int bin097 = axis2->FindBin(0.097);
  int bin112 = axis2->FindBin(0.112);
  int bin162 = axis2->FindBin(0.162);
  int bin177 = axis2->FindBin(0.177);
  int bin227 = axis2->FindBin(0.227);

  for(int ipt=0; ipt<30; ipt++)
  {
    axis1->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(2);

    double nphoton = h_1photon->GetBinContent(ipt+1);
    double npair = 0.;
    double nbgside = 0.;

    for(int ib=bin047; ib<bin097; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    for(int ib=bin177; ib<bin227; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    nbgside /= 2.;

    for(int ib=bin112; ib<bin162; ib++)
      npair += h_inv_mass->GetBinContent(ib);

    //double nsub = npair - nbgside;
    //double npion = nsub * ( 1. + MissR[ipt] );
    //double ndecay = npion * ( 1. + BR[ipt] );
    //double ndirect = nphoton - ndecay;

    yield[0*30+ipt] = nphoton;
    yield[1*30+ipt] = npair;
    yield[2*30+ipt] = nbgside;

    delete h_inv_mass;
  }

  delete h_1photon;
  return;
}

void generate_ALL(int Process = 0)
{
  gSystem->Load("libDirectPhotonPP.so");

  const int npT = 30;
  const float pTbins[npT+1] = { 0.0,
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
    12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0};

  ifstream fin("/phenix/plhf/zji/taxi/Run13pp510ERT/runlist.txt");
  int nrun = 0;
  int runnumber[1024];
  while(fin >> runnumber[nrun]) nrun++;
  fin.close();

  const int nP = 50;
  const int firstRun = Process * nP;
  const int lastRun = (Process + 1) * nP;

  // ALL[type][part][crossing][pT][run]
  // type: 0-photon, 1-pair, 2-bgside
  // part: 0-PbScW, 1-PbScE, 2-PbGlE
  // crossing: 0-even, 1-odd
  float ALL[3][3][2][npT][1024] = {};
  float eALL[3][3][2][npT][1024] = {};

  for(int ir=firstRun; ir<lastRun; ir++)
  {
    TFile *f1 = new TFile(Form("/phenix/plhf/zji/taxi/Run13pp510ERT/8847/data/DirectPhotonPP_PhotonNode-%d.root", runnumber[ir]));
    TFile *f2 = new TFile(Form("/phenix/plhf/zji/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/PhotonNode-ERT-%d.root", runnumber[ir]));
    TTree *T1 = (TTree*)f1->Get("T1");

    SpinPattern *spinpat = new SpinPattern;
    T1->SetBranchAddress("RUN/SpinPattern", &spinpat);
    T1->GetEntry(0);

    // 0-even_crossing, 1-odd_crossing
    Long64_t bbc_same[2] = {};
    Long64_t bbc_opposite[2] = {};
    for(int ib=0; ib<120; ib++)
    {
      int pattern = spinpat->get_spinpattern_blue(ib) * spinpat->get_spinpattern_yellow(ib);
      if(pattern == 1)
        bbc_same[ib%2] += spinpat->get_bbc_wide(ib);
      else if(pattern == -1)
        bbc_opposite[ib%2] += spinpat->get_bbc_wide(ib);
    }

    float pb = spinpat->get_pb();
    float py = spinpat->get_py();
    float epb = spinpat->get_pbstat();
    float epy = spinpat->get_pystat();

    if( bbc_same[0] <= 0 || bbc_opposite[0] <= 0 ||
        bbc_same[1] <= 0 || bbc_opposite[1] <= 0 ||
        pb <= 0. || py <= 0. )
      continue;

    // yield[type][pattern][part][ipt]
    // type: 0-photon, 1-pair, 2-bgside
    // pattern: 0-same, 1-opposite
    // part: 0-PbScW, 1-PbScE, 2-PbGlE
    Long64_t yield[3][2][3][npT] = {};

    for(int itype=0; itype<3; itype++)
      for(int ipat=0; ipat<2; ipat++)
        for(int ipart=0; ipart<3; ipart++)
        {
          int pattern = ( ipat==0 ? 1 : -1 );
          GetRawAsymm(f2, itype, pattern, (Long64_t*)yield[itype][ipat]);
        }

    for(int itype=0; itype<3; itype++)
      for(int ipart=0; ipart<3; ipart++)
        for(int icr=0; icr<2; icr++)
          for(int ipt=0; ipt<npT; ipt++)
          {
            float yield_same = yield[itype][0][ipart][ipt];
            float yield_opposite = yield[itype][1][ipart][ipt];
            if( yield_same <= 0. || yield_opposite <= 0. ) continue;
            float R = (float)bbc_same[icr] / bbc_opposite[icr];
            float eR = R * sqrt( 1./yield_same + 1./yield_opposite );
            ALL[itype][ipart][icr][ipt][ir] = ( yield_same - R * yield_opposite ) / ( yield_same + R * yield_opposite ) / ( pb * py );
            eALL[itype][ipart][icr][ipt][ir] = sqrt( pow( 2. * R * yield_same * yield_opposite / (pb * py) / pow(yield_same + R * yield_opposite, 2.) , 2. )
                * pow( 1./yield_same + 1./yield_opposite + (eR*eR)/(R*R), 2.)
                + ( (epb*epb)/(pb*pb) + (epy*epy)/(py*py) ) * pow(ALL[itype][ipart][icr][ipt][ir], 2.) );
          }

    delete spinpat;
    delete f1;
    delete f2;
  }

  for(int itype=0; itype<3; itype++)
    for(int ipart=0; ipart<3; ipart++)
      for(int icr=0; icr<2; icr++)
        for(int ipt=0; ipt<npT; ipt++)
        {
          ofstream fout(Form("ALL/Process%d-type%d-part%d-crossing%d-pT%0d.txt", Process, itype, ipart, icr, ipt), ofstream::trunc);
          for(int ir=firstRun; ir<lastRun; ir++)
            if( eALL[itype][ipart][icr][ipt][ir] > 0. )
              fout << runnumber[ir] << "\t" << ALL[itype][ipart][icr][ipt][ir]
                << "\t" << eALL[itype][ipart][icr][ipt][ir] << endl;
          fout.close();
        }
}
