#ifndef __PTWEIGHTS_H__
#define __PTWEIGHTS_H__

class TF1;
class TFile;
class TH2;

class PtWeights
{
  public:
    PtWeights();
    virtual ~PtWeights();

    double EvalPi0(double pt);
    double EvalPhoton(double pt);
    double Integral(double pt1, double pt2, const char *option);

  protected:
    void ReadWeights();

    TF1 *cross_pi0;
    TF1 *cross_ph;
    TFile *f_mb;
    TH2 *h2_mb;
};

#endif /* __PTWEIGHTS_H__ */
