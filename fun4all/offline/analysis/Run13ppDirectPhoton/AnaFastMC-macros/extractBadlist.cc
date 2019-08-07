#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
  if(argc < 3)
  {
    cout << "Usage: " << argv[0] << " <goodlist.txt> <nevents>" << endl;
    return 1;
  }

  int nproc;
  stringstream ss;
  ss << argv[2];
  ss >> nproc;

  ifstream fin(argv[1]);
  if(!fin.is_open() || nproc < 0)
  {
    cout << "Wrong input" << endl;
    return 1;
  }

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

  fin.close();
  return 0;
}
