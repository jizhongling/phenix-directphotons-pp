void extractBadlist()
{
  const int nproc = 12000;
  ifstream fin("/phenix/spin/phnxsp01/zji/data/pythiaToHisto/aa_goodlist.txt");
  vector<int> goodlist(nproc);
  int proc;
  while(fin >> proc)
    goodlist.push_back(proc);

  for(int proc=0; proc<nproc; proc++)
  {
    bool isbad = true;
    for(vector<int>::iterator it=goodlist.begin(); it != goodlist.end(); it++)
      if(proc == *it)
      {
        isbad = false;
        break;
      }
    if(isbad)
      cout << proc << " ";
  }
  cout << endl;
}
