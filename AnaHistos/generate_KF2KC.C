void generate_KF2KC()
{
  const Int_t N = 1e6;
  Int_t vkf[N];
  Int_t vkc[N];
  Int_t idx = 0;
  TPythia6 *tpythia6 = new TPythia6();

  for(Int_t kf=1; kf<N; kf++)
  {
    if(kf > 1000 && kf/10%10 == 0)
      continue;
    Int_t kc = tpythia6->Pycomp(kf);
    if(kc != 0)
    {
      bool has = false;
      for(Int_t i=0; i<idx; i++)
        if(kc == vkc[i])
          has = true;
      if(!has)
      {
        vkf[idx] = kf;
        vkc[idx] = kc;
        idx++;
      }
    }
  }

  for(Int_t i=0; i<idx; i++)
    cout << "mdcy " << vkc[i] << " 1 0\t// " << vkf[i] << endl;
}
