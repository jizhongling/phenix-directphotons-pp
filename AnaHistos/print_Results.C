#include "GlobalVars.h"
#include "QueryTree.h"
#include "IsoPhotonALL.h"

void print(const char *value_name, const char *sys_name, int value_part, int sys_part, bool ispol)
{
  QueryTree *qt_value = new QueryTree(Form("data/%s.root",value_name));
  QueryTree *qt_sys = new QueryTree(Form("data/%s.root",sys_name));

  for(int ipt=0; ipt<(ispol?npT_pol:npT); ipt++)
  {
    double xpt, value, stat, sys;
    qt_value->Query(ipt, value_part, xpt, value, stat);
    qt_sys->Query(ipt, sys_part, xpt, value, sys);
    if(ispol)
      printf("%.1f-%.1f", pTbin_pol[ipt], pTbin_pol[ipt+1]);
    else
      printf("%.1f-%.1f", pTbin[ipt], pTbin[ipt+1]);
    cout << " & " << value << " & " << stat << " & " << sys << " \\\\" << endl;
  }
  cout << value_name << (ispol?" beam ":" iso ") << sys_part << endl << endl;

  return;
}

void print_Results()
{
  for(int iso=0; iso<2; iso++)
  {
    char *name = iso ? "CrossSection-isophoton" : "CrossSection-photon";
    print(name, "CrossSection-syserr", 3, iso, false);
  }
  for(int beam=0; beam<3; beam++)
  {
    int igr = beam + ngr_photon*2;
    print("IsoPhotonALL", "IsoPhotonALL-syserr", igr, beam, true);
  }
}
