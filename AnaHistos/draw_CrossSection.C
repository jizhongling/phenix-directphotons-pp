void GenerateGraph(TFile *f, Int_t part, Int_t ispion, Int_t gn, Double_t **gx, Double_t **gy, Double_t **egy)
{
  const Double_t PI = 3.1415927;
  const Double_t BBCCount = 4.193255 * pow(10,12);  // 4193255288313
  const Double_t BBCRatio = 0.54;
  const Double_t eBBCRatio = 0.09;
  const Double_t BBCCross = 32.5 * pow(10,9);
  const Double_t eBBCCross = 3.25 * pow(10,9);
  const Double_t BR[30] = {0.0444823,0.155126,0.189877,0.209775,0.227416,0.235346,0.242619,0.258304,0.254693,0.219828,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313,0.243313};
  const Double_t eBR[30] = {1.60171e-05,0.00010065,0.000293243,0.000651539,0.00129034,0.00232753,0.00398618,0.00663705,0.0101993,0.0144502,0.0228892,0.0291354,0.0348798,0.0512812,0.0699279,0.078914,0.111523,0.142477,0.247947,0.213559,0.156142,0.154185,0.30816,0.8415,0.8415};

  if(ispion == 0)
  {
    if(part == 0)
    {
      const Int_t sector_low = 1;
      const Int_t sector_high = 4;
      const Double_t MissR[30] = {1.2553, 1.17067, 1.59935, 1.42668, 1.25836, 1.33632, 1.23448, 1.24616, 1.21511, 1.20831, 1.06156, 1.11788, 1.03933, 1.09151, 0.98973, 0.823309, 0.869805, 0.935766, 0.840094, 0.731882, 0.751506, 0.723598, 0.707506, 0.828866, 0.840945, 1.15116, 1.56706, 2.31187, 3.97172, 6.62838};
      const Double_t eMissR[30] = {0.35373, 0.0340873, 0.0511709, 0.0507457, 0.0502399, 0.056308, 0.0576832, 0.0617297, 0.0631708, 0.064343, 0.060742, 0.0647261, 0.0622611, 0.0667717, 0.0612381, 0.0529919, 0.0561243, 0.0625699, 0.0559797, 0.0509192, 0.0283824, 0.0305927, 0.0322388, 0.0403963, 0.0435559, 0.0588781, 0.0809105, 0.128253, 0.238921, 0.571009};
      const Double_t Acceptance[30] = {0,0,0.0195455,0.0198703,0.0176629,0.0171196,0.0144225,0.0129019,0.0128993,0.0128214,0.011992,0.0127937,0.0124924,0.0125693,0.0133918,0.0128314,0.0142356,0.0133699,0.0142002,0.0142407,0.0163306,0.0171801,0.0180004,0.0189751,0.0209591,0.0272184,0.0381716,0.0534764,0.111023,0.271603};
      const Double_t eAcceptance[30] = {0,0,6.0226e-05,6.47864e-05,6.04221e-05,6.08869e-05,5.33453e-05,4.90899e-05,5.07674e-05,5.19692e-05,5.00023e-05,5.46465e-05,5.50648e-05,5.6422e-05,6.21171e-05,6.07512e-05,6.87508e-05,6.61692e-05,7.18825e-05,7.37803e-05,4.47378e-05,5.14826e-05,5.93346e-05,6.94026e-05,8.51554e-05,0.000125162,0.000203089,0.000342292,0.000933821,0.00400457};
      const Double_t TrigE[30] = {8.09532e-06,1.88351e-05,9.80708e-05,0.000389597,0.00105885,0.00353447,0.0152754,0.073618,0.227824,0.470076,0.679126,0.884265,0.909871,0.953846,0.952096,0.954128,0.897059,0.958333,0.965517,1,0.953488,0.916667,1,1,1,1,1,1,1,1};
      const Double_t eTrigE[30] = {2.47946e-06,7.6374e-07,5.0765e-06,2.25507e-05,7.22226e-05,0.000230904,0.000770765,0.00250532,0.00587997,0.00991908,0.012883,0.0127204,0.0153001,0.016987,0.0227901,0.0298569,0.0510686,0.0523567,0.074913,0.0615421,0.0581073,0.166603,0.264349,0.369032,0.417333,0.601879,0.8415,0.8415,0.601879,0.8415};
    }
    else if(part == 1)
    {
      const Int_t sector_low = 5;
      const Int_t sector_high = 6;
      const Double_t MissR[30] = {0.978365, 1.24938, 1.8402, 1.56844, 1.1526, 1.26383, 1.3897, 1.16472, 1.26581, 1.05511, 1.08663, 0.955101, 1.10811, 1.08163, 0.95795, 0.855263, 0.732881, 0.790301, 0.704059, 0.742194, 0.678615, 0.677137, 0.660269, 0.75477, 0.854586, 1.29879, 1.72092, 2.77472, 4.19073, 8.35532};
      const Double_t eMissR[30] = {0.370777, 0.057034, 0.0872807, 0.0795781, 0.0653519, 0.0808292, 0.0906953, 0.0822108, 0.0943503, 0.0855029, 0.0872746, 0.0829433, 0.0915098, 0.0910154, 0.0844898, 0.0781328, 0.0706467, 0.0806208, 0.0685976, 0.0743618, 0.0380066, 0.0419326, 0.0433945, 0.0544783, 0.0642592, 0.101866, 0.132278, 0.229672, 0.37979, 1.08071};
      const Double_t Acceptance[30] = {0,0,0.00935206,0.010297,0.00907584,0.00749923,0.00784875,0.00692969,0.00646637,0.00601953,0.00624646,0.00583564,0.00689024,0.00673594,0.0068271,0.00684924,0.00675814,0.00602485,0.00721831,0.00678838,0.00780776,0.00859304,0.0090514,0.00893903,0.0102763,0.0125893,0.0188214,0.0272568,0.0554809,0.164154};
      const Double_t eAcceptance[30] = {0,0,2.88167e-05,3.3573e-05,3.10471e-05,2.66714e-05,2.90306e-05,2.63665e-05,2.54495e-05,2.4399e-05,2.60454e-05,2.49261e-05,3.03712e-05,3.02367e-05,3.16671e-05,3.24281e-05,3.26385e-05,2.98176e-05,3.65396e-05,3.51701e-05,2.13893e-05,2.57503e-05,2.98361e-05,3.26951e-05,4.17519e-05,5.78912e-05,0.000100138,0.000174465,0.000466654,0.00242031};
      const Double_t TrigE[30] = {8.45417e-06,2.15626e-05,0.00012377,0.00043417,0.00106439,0.0057713,0.0278967,0.0999358,0.240773,0.414296,0.604027,0.805699,0.828,0.888889,0.95122,0.941176,0.934783,0.96875,0.944444,0.909091,0.9375,0.888889,1,1,1,1,1,1,1,1};
      const Double_t eTrigE[30] = {3.86319e-06,1.17435e-06,7.86683e-06,3.27608e-05,0.000100533,0.000399679,0.00139129,0.00394678,0.00820872,0.0137363,0.0187326,0.0222699,0.0272703,0.0315277,0.0369117,0.0440953,0.0594246,0.068259,0.116475,0.179384,0.076593,0.21166,0.601879,0.601879,0.8415,0.8415,0.8415,0.8415,0.8415,0.8415};
    }
    else if(part == 2)
    {
      const Int_t sector_low = 7;
      const Int_t sector_high = 8;
      const Double_t MissR[30] = {1.61424, 1.38534, 1.80124, 1.5525, 1.34054, 1.28686, 1.23901, 1.14929, 1.15866, 1.23586, 0.974419, 1.02197, 0.917833, 0.795642, 0.789946, 0.84228, 0.789443, 0.845898, 0.707791, 0.591573, 0.589901, 0.56227, 0.57542, 0.573242, 0.508064, 0.554727, 0.681324, 1.1311, 1.6023, 5.19963};
      const Double_t eMissR[30] = {0.664857, 0.0581373, 0.0791804, 0.0759765, 0.070702, 0.0743297, 0.0787384, 0.0753755, 0.0807072, 0.0879502, 0.0779303, 0.0809514, 0.073861, 0.0676877, 0.0646398, 0.0698503, 0.0668627, 0.071638, 0.0614171, 0.0522186, 0.0289569, 0.0290219, 0.0308298, 0.033052, 0.0317432, 0.0368717, 0.04619, 0.0830421, 0.118913, 0.774534};
      const Double_t Acceptance[30] = {0,0,0.0103337,0.0103094,0.00966938,0.00860223,0.00748387,0.00698105,0.00643502,0.00673056,0.0055803,0.00591335,0.00605819,0.00564972,0.00693004,0.0068375,0.00687547,0.00733694,0.00721827,0.00809117,0.00889558,0.00994597,0.0112441,0.0115621,0.0136592,0.0147227,0.0196718,0.0214841,0.0420287,0.0727782};
      const Double_t eAcceptance[30] = {0,0,3.18416e-05,3.36135e-05,3.30775e-05,3.05943e-05,2.7681e-05,2.6562e-05,2.53261e-05,2.7281e-05,2.32678e-05,2.52581e-05,2.67036e-05,2.53608e-05,3.21446e-05,3.23725e-05,3.32052e-05,3.63113e-05,3.65394e-05,4.19198e-05,2.43694e-05,2.98045e-05,3.70637e-05,4.22891e-05,5.54964e-05,6.77014e-05,0.000104662,0.000137515,0.000353507,0.00107306};
      const Double_t TrigE[30] = {3.86368e-05,6.81334e-05,0.000125886,0.000260003,0.000887477,0.0019329,0.0187945,0.0556734,0.128918,0.239865,0.351272,0.488696,0.662824,0.638393,0.690141,0.681818,0.760563,0.790698,0.709677,0.757576,0.808511,0.866667,0.625,0.5,0.5,1,1,1,1,1};
      const Double_t eTrigE[30] = {7.23467e-06,2.04833e-06,8.44487e-06,2.56335e-05,8.06052e-05,0.000140612,0.00106058,0.00276853,0.00577574,0.0105686,0.0155742,0.0217299,0.027293,0.0349367,0.0436609,0.0573316,0.06149,0.0808377,0.103954,0.0976443,0.0751301,0.149614,0.235052,0.417333,0.417333,0.601879,0.8415,0.8415,0.8415,0.8415};
    }
    else
    {
      cerr << "Wrong part " << part << endl;
      exit(1);
    }
  }
  else if(ispion == 1)
  {
    if(part == 0)
    {
      const Int_t sector_low = 1;
      const Int_t sector_high = 4;
      const Double_t MissR[30] = {0, 0, 8.05882, 2.79473, 2.11713, 1.52164, 1.41379, 1.16843, 1.14565, 1.07649, 1.05306, 1.01185, 0.871967, 0.847773, 0.817959, 0.784861, 0.779126, 0.75117, 0.745971, 0.737539, 0.661434, 0.635864, 0.581672, 0.644476, 0.636364, 0.751874, 0.940488, 1.06953, 1.39528, 1.63534};
      const Double_t eMissR[30] = {0, 0, 0.598214, 0.141323, 0.0982971, 0.0669906, 0.0606414, 0.0493818, 0.0490194, 0.0456639, 0.0444753, 0.0415176, 0.0369586, 0.0356148, 0.0348409, 0.03341, 0.0334887, 0.0320324, 0.031616, 0.031592, 0.0144172, 0.0145506, 0.014518, 0.0174388, 0.0187944, 0.0227948, 0.0304678, 0.0346741, 0.0455753, 0.053406};
      const Double_t Acceptance[30] = {0,0,0.0146609,0.0390878,0.0462119,0.0534222,0.0565653,0.0584277,0.0581461,0.0609997,0.0638506,0.0666351,0.0681279,0.0679608,0.0683042,0.0658338,0.0715757,0.070512,0.0725108,0.0726908,0.0723308,0.0738418,0.0610559,0.0483832,0.0381794,0.0301752,0.023007,0.017892,0.0148446,0.00933417};
      const Double_t eAcceptance[30] = {0,0,0.000868189,0.00151891,0.00163354,0.00179564,0.00183798,0.00189551,0.00191905,0.0020071,0.00205222,0.00208959,0.00213239,0.00215451,0.00213759,0.00209551,0.00225711,0.00216955,0.00225102,0.00222605,0.00112079,0.00124284,0.00108269,0.000938835,0.000825074,0.000720094,0.000617,0.000518839,0.00052161,0.000413444};
      const Double_t TrigE[30] = {0,3.14701e-05,7.41417e-05,0.000262871,0.000382521,0.00138882,0.00290081,0.00895857,0.0362069,0.0880551,0.179085,0.291855,0.417004,0.589474,0.707317,0.855556,0.810345,0.847826,0.88,0.885246,1,1,1,1,1,1,1,1,1,1};
      const Double_t eTrigE[30] = {0.00052301,1.13783e-05,2.1221e-05,5.96839e-05,0.000117131,0.000321933,0.000686935,0.00164447,0.00430998,0.0085611,0.0149588,0.0232406,0.0337037,0.0387411,0.046783,0.0470126,0.0658858,0.0723547,0.10317,0.0563623,0.185081,0.458818,0.601879,0.601879,0.8415,0.8415,0.8415,0.8415,0.8415,0.8415};
    }
    else if(part == 1)
    {
      const Int_t sector_low = 5;
      const Int_t sector_high = 6;
      const Double_t MissR[30] = {0, 0, 7.32075, 3.36325, 2.04491, 1.67102, 1.46506, 1.24307, 1.11321, 0.938959, 0.928016, 0.941767, 0.940476, 0.84029, 0.864504, 0.762478, 0.836996, 0.714984, 0.723708, 0.727928, 0.616941, 0.565986, 0.532725, 0.581805, 0.669202, 0.838493, 1.12797, 1.74442, 2.64952, 4.32618};
      const Double_t eMissR[30] = {0, 0, 0.758065, 0.250424, 0.136537, 0.107952, 0.0932862, 0.0771051, 0.0702263, 0.0571715, 0.0589999, 0.0605976, 0.0601745, 0.0529764, 0.0554625, 0.0480936, 0.0530664, 0.0446883, 0.0471555, 0.0476059, 0.0205031, 0.020049, 0.0209296, 0.0254758, 0.0325855, 0.0455498, 0.0662427, 0.108993, 0.176328, 0.314472};
      const Double_t Acceptance[30] = {0,0,0.00674877,0.0155559,0.0208567,0.0258659,0.0263526,0.027224,0.0300626,0.0337466,0.0323423,0.0343707,0.0304741,0.0348289,0.0364121,0.0339585,0.0400664,0.0372991,0.0371205,0.0335042,0.0373322,0.0380132,0.0332746,0.0248774,0.0186953,0.0130948,0.00886884,0.00604235,0.00511665,0.00368615};
      const Double_t eAcceptance[30] = {0,0,0.000624577,0.000991379,0.00114976,0.00129835,0.00129049,0.00131802,0.00141121,0.00151739,0.00151651,0.00155614,0.00145536,0.00158208,0.0016194,0.00153756,0.0017122,0.00159847,0.00163698,0.00152833,0.000823994,0.000919641,0.000836123,0.000696315,0.000612222,0.000501825,0.000404434,0.000324311,0.000331413,0.000278556};
      const Double_t TrigE[30] = {0,1.7545e-05,7.37063e-05,0.000278379,0.000436904,0.00129429,0.00430515,0.0151783,0.0435967,0.106529,0.188017,0.358885,0.354651,0.526718,0.682927,0.829787,0.756098,0.857143,0.904762,0.84375,0.888889,1,1,1,1,1,1,1,1,1};
      const Double_t eTrigE[30] = {0.000869721,1.38804e-05,2.96093e-05,8.23269e-05,0.000166066,0.000411081,0.00104254,0.00266317,0.0060084,0.0114964,0.0194737,0.0305431,0.0402014,0.0474924,0.0596374,0.0732132,0.0856243,0.119782,0.112122,0.0919399,0.21166,0.30816,0.8415,0.8415,0.8415,0.8415,0.8415,0.8415,0.8415,0.8415};
    }
    else if(part == 2)
    {
      const Int_t sector_low = 7;
      const Int_t sector_high = 8;
      const Double_t MissR[30] = {0, 0, 10.8, 3.42339, 2.34049, 1.73031, 1.51584, 1.21373, 1.20849, 1.10761, 0.866352, 0.803954, 0.894198, 0.746106, 0.710937, 0.740171, 0.623145, 0.664168, 0.656394, 0.615504, 0.579342, 0.539318, 0.492479, 0.434566, 0.45366, 0.427498, 0.441251, 0.536877, 0.651888, 0.801914};
      const Double_t eMissR[30] = {0, 0, 1.26214, 0.247104, 0.154864, 0.106184, 0.0928873, 0.0725834, 0.0717803, 0.0658102, 0.0504215, 0.0488803, 0.0537627, 0.0450472, 0.0435956, 0.0469228, 0.0387386, 0.0407075, 0.04093, 0.0392635, 0.0180835, 0.0172653, 0.0162252, 0.0148895, 0.0160251, 0.0164251, 0.0179082, 0.0234147, 0.029734, 0.0371854};
      const Double_t Acceptance[30] = {0,0,0.00495649,0.0144047,0.0192193,0.0249476,0.0264911,0.0275352,0.0300894,0.0309558,0.0335992,0.0337508,0.0316126,0.0348021,0.0335877,0.0336534,0.0388197,0.0343989,0.0369946,0.0307659,0.038342,0.0402307,0.0394662,0.0375388,0.0357642,0.0303404,0.0259229,0.0187823,0.0154237,0.00975407};
      const Double_t eAcceptance[30] = {0,0,0.000524366,0.00092323,0.00106765,0.00122427,0.00124626,0.00126861,0.00133072,0.00134878,0.00140843,0.00141235,0.00134647,0.00144594,0.00141343,0.00140528,0.00152915,0.00140927,0.00148437,0.0012992,0.000769355,0.000816018,0.000789077,0.000763113,0.00074863,0.00068482,0.000629957,0.000527182,0.000491757,0.000395857};
      const Double_t TrigE[30] = {0,0.000113637,0.000166582,0.000183883,0.000315468,0.0011467,0.00392341,0.00687799,0.0236524,0.0510397,0.0708402,0.144357,0.196,0.257143,0.336538,0.426471,0.444444,0.606061,0.619048,0.679245,0.6875,0.833333,0.333333,0.666667,1,1,1,1,1,1};
      const Double_t eTrigE[30] = {0.00114487,2.63621e-05,4.04438e-05,6.99055e-05,0.000144123,0.000379018,0.000950324,0.00174744,0.00412902,0.0077023,0.0120125,0.0203188,0.0283943,0.0373284,0.0526449,0.0680235,0.0772817,0.102252,0.13297,0.0766927,0.156142,0.1791,0.414672,0.414672,0.8415,0.8415,0.8415,0.8415,0.8415,0.8415};
    }
    else
    {
      cerr << "Wrong part " << part << endl;
      exit(1);
    }
  }
  else
  {
    cerr << "Wrong ispion " << ispion << endl;
    exit(1);
  }

  THnSparse *hn_1photon = (THnSparse*)f->Get("hn_1photon");
  hn_1photon->GetAxis(0)->SetRange(sector_low, sector_high);
  TH1 *h_1photon = hn_1photon->Projection(1);

  TH2 *h2_sig_extra = (TH2*)f->Get("h2_sig_extra");
  TH1 *h_sig_extra = h2_sig_extra->ProjectionY("h_sig_extra", sector_low, sector_high);

  TH2 *h2_bg_extra = (TH2*)f->Get("h2_bg_extra");
  TH1 *h_bg_extra = h2_bg_extra->ProjectionY("h_bg_extra", sector_low, sector_high);

  THnSparse *hn_2photon = (THnSparse*)f->Get("hn_2photon");

  TAxis *axis0 = hn_2photon->GetAxis(0);
  if(ispion == 0)
    TAxis *axis1_2 = hn_2photon->GetAxis(1);
  else if(ispion == 1)
    TAxis *axis1_2 = hn_2photon->GetAxis(2);
  TAxis *axis3 = hn_2photon->GetAxis(3);

  axis0->SetRange(sector_low, sector_high);

  Int_t bin047 = axis3->FindBin(0.047);
  Int_t bin067 = axis3->FindBin(0.067);
  Int_t bin087 = axis3->FindBin(0.087);
  Int_t bin097 = axis3->FindBin(0.097);
  Int_t bin112 = axis3->FindBin(0.112);
  Int_t bin162 = axis3->FindBin(0.162);
  Int_t bin177 = axis3->FindBin(0.177);
  Int_t bin187 = axis3->FindBin(0.187);
  Int_t bin212 = axis3->FindBin(0.212);
  Int_t bin227 = axis3->FindBin(0.227);

  const Int_t nData = 256;
  vector<Double_t> x(nData), y(nData), sigma_y(nData);

  TCanvas *c = new TCanvas("c", "Canvas", 2400, 2000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c->Divide(6,5);

  Int_t ipad = 1;
  for(Int_t ipt=0; ipt<gn; ipt++)
  {
    c->cd(ipad++);

    char buf[100];
    Double_t low = axis1_2->GetBinLowEdge(ipt+1);
    Double_t high = axis1_2->GetBinLowEdge(ipt+2);
    sprintf(buf, "p_{T} %3.1f-%3.1f", low, high);

    axis1_2->SetRange(ipt+1,ipt+1);
    TH1 *h_inv_mass = hn_2photon->Projection(3);
    h_inv_mass->SetTitle(buf);

    Double_t nphoton = h_1photon->GetBinContent(ipt+1);
    Double_t nsig_extra = h_sig_extra->GetBinContent(ipt+1);
    Double_t nbg_extra = h_bg_extra->GetBinContent(ipt+1);
    Double_t nextra = nsig_extra - nbg_extra/2.;
    Double_t npair = 0.;
    Double_t nbgfit = 0.;
    Double_t nbgside = 0.;
    Double_t nbggpr = 0.;
    Double_t dnbggpr = 0.;

    x.clear();
    y.clear();
    sigma_y.clear();

    for(Int_t ib=bin067; ib<bin087; ib++)
    {
      Double_t xx = axis3->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    for(Int_t ib=bin187; ib<bin212; ib++)
    {
      Double_t xx = axis3->GetBinCenter(ib);
      Double_t yy = h_inv_mass->GetBinContent(ib);
      Double_t sigma_yy = h_inv_mass->GetBinError(ib);
      x.push_back(xx);
      y.push_back(yy);
      sigma_y.push_back(sigma_yy);
    }

    BgGPR(x, y, sigma_y, nbggpr, dnbggpr);
    nbggpr /= 0.001;
    dnbggpr /= 0.001;

    TF1 *fn1 = new TF1("fn1", "gaus", 0., 0.5);
    TF1 *fn2 = new TF1("fn2", "gaus(0)+pol3(3)", 0., 0.5);
    TF1 *fn3 = new TF1("fn3", "pol3", 0., 0.5);

    Double_t par[10];
    h_inv_mass->Fit(fn1, "Q0", "", 0.112, 0.162);
    fn2->SetParameters( fn1->GetParameters() );
    h_inv_mass->Fit(fn2, "Q", "", 0.047, 0.227);
    fn2->GetParameters(par);
    fn3->SetParameters(par[3], par[4], par[5], par[6]);

    for(Int_t ib=bin047; ib<bin097; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    for(Int_t ib=bin177; ib<bin227; ib++)
      nbgside += h_inv_mass->GetBinContent(ib);
    nbgside /= 2.;

    for(Int_t ib=bin112; ib<bin162; ib++)
    {
      npair += h_inv_mass->GetBinContent(ib);
      Double_t bincenter = axis3->GetBinCenter(ib);
      nbgfit += fn3->Eval(bincenter);
    }

    Double_t nfit = npair - nbgfit;
    Double_t enfit = sqrt(1./npair + 1./nbgfit);
    Double_t nsub = npair - nbgside;
    Double_t ndiff = fabs(nbgfit - nbgside - nextra);
    Double_t ngpr = npair - nbggpr;
    Double_t nerror = fabs(nfit - ngpr);

    if(nfit > 0. && nfit < npair)
      Double_t npion = nfit * ( 1. + MissR[ipt] );
    else
      Double_t npion = nsub * ( 1. + MissR[ipt] );
    Double_t enpion = npion * sqrt( pow(enfit/nfit,2.) + pow(eMissR[ipt]/(1.+MissR[ipt]),2.) );
    Double_t ndecay = npion * ( 1. + BR[ipt] );
    Double_t endecay = ndecay * sqrt( pow(enpion/npion,2.) + pow(eBR[ipt]/(1.+BR[ipt]),2.) );
    Double_t ndirect = nphoton - ndecay;
    Double_t endirect = sqrt( 1./nphoton + pow(endecay/ndecay,2.) );

    if(ispion == 0)
    {
      gx[0][ipt] = (low + high) / 2.;
      gy[0][ipt] = ndirect / nphoton;
      egy[0][ipt] = gy[0][ipt] * sqrt( pow(endirect/ndirect,2.) + 1./nphoton );

      gx[1][ipt] = (low + high) / 2.;
      gy[1][ipt] = BBCCross * (ndirect/(BBCCount/BBCRatio)) / (2*PI*gx[1][ipt]) / ((high-low)*0.8) / Acceptance[ipt] / TrigE[ipt];
      egy[1][ipt] = gy[1][ipt] * sqrt( pow(eBBCCross/BBCCross,2.) + pow(endirect/ndirect,2.) + 1./BBCCount + pow(eBBCRatio/BBCRatio,2.)
          + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.) );

      cout << "pT=" << low << "-" << high << "\tnphoton=" << nphoton
        << "\tnpair=" << npair << "\tnfit=" << nfit << "\tnsub=" << nsub
        << "\tndiff=" << ndiff/nfit*100. << "%\tngpr=" << ngpr
        << "\tnerror=" << nerror/nfit*100. << "%" << endl;
    }
    else if(ispion == 1)
    {
      gx[2][ipt] = (low + high) / 2.;
      gy[2][ipt] = BBCCross * (npion/2./(BBCCount/BBCRatio)) / (2*PI*gx[2][ipt]) / ((high-low)*0.8) / Acceptance[ipt] / TrigE[ipt];
      egy[2][ipt] = gy[2][ipt] * sqrt( pow(eBBCCross/BBCCross,2.) + pow(enpion/npion,2.) + 1./BBCCount + pow(eBBCRatio/BBCRatio,2.)
          + pow(eAcceptance[ipt]/Acceptance[ipt],2.) + pow(eTrigE[ipt]/TrigE[ipt],2.) );
    }
  }


  if(ispion == 0)
  {
    if(part == 0)
    {
      c->Print("CrossSection-PbScW.pdf");
      delete c;
    }
    if(part == 1)
    {
      c->Print("CrossSection-PbScE.pdf");
      delete c;
    }
    else if(part == 2)
    {
      c->Print("CrossSection-PbGlE.pdf");
      delete c;
    }
  }
  else if(ispion == 1)
  {
    if(part == 0)
    {
      c->Print("CrossSection-PbScW-pion.pdf");
      delete c;
    }
    if(part == 1)
    {
      c->Print("CrossSection-PbScE-pion.pdf");
      delete c;
    }
    else if(part == 2)
    {
      c->Print("CrossSection-PbGlE-pion.pdf");
      delete c;
    }
  }

  return;
}

TGraphErrors **CreateGraph(TFile *f, Int_t part)
{
  const Int_t gn = 30;
  Double_t *gx[3];
  Double_t *gy[3];
  Double_t *egy[3];
  for(Int_t i=0; i<3; i++)
  {
    gx[i] = new Double_t[gn];
    gy[i] = new Double_t[gn];
    egy[i] = new Double_t[gn];
    for(Int_t j=0; j<gn; j++)
    {
      gx[i][j] = 0.;
      gy[i][j] = 0.;
      egy[i][j] = 0.;
    }
  }

  for(Int_t ispion=0; ispion<2; ispion++)
    GenerateGraph(f, part, ispion, gn, gx, gy, egy);

  TGraphErrors **graph = new TGraphErrors*[3];
  for(Int_t i=0; i<3; i++)
    graph[i] = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
  return graph;
}

void draw_CrossSection()
{
  gSystem->Load("libGausProc.so");
  gROOT->ProcessLine(".L BgGPR.C");
  gROOT->ProcessLine(".L Chi2Fit.C");

  TFile *f = new TFile("/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/PhotonNode-macros/histos/total.root");

  TCanvas *c0 = new TCanvas("c0", "Canvas", 1800, 1200);
  gStyle->SetOptStat(0);
  c0->Divide(3,2);

  TGraphErrors **gr[3];
  for(Int_t part=0; part<3; part++)
  {
    cout << "part " << part << endl;
    gr[part] = CreateGraph(f, part);
    gr[part][0]->SetTitle("DirectPhoton fraction");
    gr[part][1]->SetTitle("DirectPhoton Cross Section");
    gr[part][2]->SetTitle("#pi^{0} Cross Section");
    for(Int_t i=0; i<3; i++)
    {
      c0->cd(i+1);
      if(i==0)
      {
        gr[part][i]->GetYaxis()->SetTitle("Ratio");
        gr[part][i]->GetYaxis()->SetRangeUser(0., 1.2);
      }
      else
      {
        gPad->SetLogy();
        gr[part][i]->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3} [pb*GeV^{-2}*c^{-3}]");
        gr[part][i]->GetYaxis()->SetRangeUser(0.01, pow(10,6));
      }
      gr[part][i]->GetXaxis()->SetTitle("p_{T} [GeV]");
      gr[part][i]->GetYaxis()->SetTitleOffset(1.2);
      gr[part][i]->GetXaxis()->SetRangeUser(0., 30.);
      gr[part][i]->SetMarkerColor(part+1);
      gr[part][i]->SetMarkerStyle(part+20);
      gr[part][i]->SetMarkerSize(1.);
      if(part == 0)
      {
        gr[part][i]->Draw("AP");
      }
      else
      {
        gr[part][i]->Draw("P");
      }
    }
  }

  for(Int_t i=0; i<3; i++)
  {
    c0->cd(i+1);
    TLegend *leg = new TLegend(0.4, 0.7, 0.7, 0.9);
    leg->AddEntry(gr[0][i], "PbScW", "P");
    leg->AddEntry(gr[1][i], "PbScE", "P");
    leg->AddEntry(gr[2][i], "PbGlE", "P");
    leg->Draw();
  }

  const Int_t gn = 30;
  Double_t *gx[3];
  Double_t *gy[3];
  Double_t *egy[3];
  for(Int_t i=0; i<3; i++)
  {
    gx[i] = new Double_t[gn];
    gy[i] = new Double_t[gn];
    egy[i] = new Double_t[gn];
    for(Int_t j=0; j<gn; j++)
    {
      gx[i][j] = 0.;
      gy[i][j] = 0.;
      egy[i][j] = 0.;
    }
  }

  for(Int_t i=0; i<3; i++)
  {
    Double_t xx[gn][3];
    Double_t yy[gn][3];
    Double_t eyy[gn][3];
    for(Int_t ipt=0; ipt<gn; ipt++)
      for(Int_t part=0; part<3; part++)
      {
        gr[part][i]->GetPoint(ipt, xx[ipt][part], yy[ipt][part]);
        eyy[ipt][part] = gr[part][i]->GetErrorY(ipt);
        gx[i][ipt] = xx[ipt][part];
      }
    for(Int_t ipt=0; ipt<gn; ipt++)
      Chi2Fit(3, (Double_t*)yy[ipt], (Double_t*)eyy[ipt], gy[i][ipt], egy[i][ipt]);
  }

  for(Int_t i=0; i<3; i++)
  {
    c0->cd(i+4);
    TGraphErrors *grt = new TGraphErrors(gn, gx[i], gy[i], 0, egy[i]);
    if(i==0)
      grt->SetTitle("DirectPhoton fraction");
    else if(i==1)
      grt->SetTitle("DirectPhoton Cross Section");
    else
      grt->SetTitle("#pi^{0} Cross Section");
    if(i==0)
    {
      grt->GetYaxis()->SetTitle("Ratio");
      grt->GetYaxis()->SetRangeUser(0., 1.2);
    }
    else
    {
      gPad->SetLogy();
      grt->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3} [pb*GeV^{-2}*c^{-3}]");
      grt->GetYaxis()->SetRangeUser(0.1, pow(10,5));
    }
    grt->GetXaxis()->SetTitle("p_{T} [GeV]");
    grt->GetYaxis()->SetTitleOffset(1.2);
    grt->GetXaxis()->SetRangeUser(0., 30.);
    grt->SetMarkerColor(1);
    grt->SetMarkerStyle(20);
    grt->SetMarkerSize(1.);
    grt->Draw("AP");
  }

  c0->Print("CrossSection.pdf");
}
