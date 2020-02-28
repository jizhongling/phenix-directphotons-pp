// r[i] = Sys[i]/Std[i]
// sum = sqrt( (r[]-1)*(r[]-1) )
void SumSysErr(int n, double r[], double er[], double &sum, double &esum)
{
  sum = 0.;
  esum = 0.;
  for(int i=0; i<n; i++)
  {
    sum += (r[i]-1)*(r[i]-1);
    esum += (r[i]-1)*(r[i]-1)*er[i]*er[i];
  }
  esum /= sum;
  sum = sqrt(sum);
  esum = sqrt(esum);

  return;
}
