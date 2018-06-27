void generate_KF2KC()
{
  const int N = 1e6;
  int vkf[N];
  int vkc[N];
  int idx = 0;
  TPythia6 *tpythia6 = new TPythia6();

  for(int kf=1; kf<N; kf++)
  {
    if(kf > 1000 && kf/10%10 == 0)
      continue;
    int kc = tpythia6->Pycomp(kf);
    if(kc != 0)
    {
      bool has = false;
      for(int i=0; i<idx; i++)
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

  for(int i=0; i<idx; i++)
    cout << "mdcy " << vkc[i] << " 1 0\t// " << vkf[i] << endl;
}
