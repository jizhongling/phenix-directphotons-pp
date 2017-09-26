from ROOT import gROOT, TCanvas, THnSparse, TF1, TFile
import matplotlib.pyplot as plt
import numpy as np

gROOT.Reset()

# open input file
f = TFile("data/DirectPhotonPP-Run13pp510ERT.root")
f.ls()

# open ROOT::TBrowser
#tb = TBrowser()

# get data histogram
hdata = f.Get("inv_mass_2photon")

# select cut level
# 1 = direct photon with isolation cut
# 2 = direct photon w/o  isolation cut
hdata.GetAxis(3).SetRange(2,2)

# scan all pT bins
# select pT bins
ptbin1 = 7
ptbin2 = 7

ptbins = range(7,24)

# create empty list to store histograms
hists_inv_mass = []

for counter, ptbin1 in enumerate(ptbins):

    # Set pT bin range for input histogram
    ptbin2 = ptbin1
    hdata.GetAxis(0).SetRange(ptbin1,ptbin2)

    # get pT range for this pT bin
    ptmin = hdata.GetAxis(0).GetBinCenter(ptbin1) - 0.5 * hdata.GetAxis(0).GetBinWidth(ptbin1)
    ptmax = hdata.GetAxis(0).GetBinCenter(ptbin2) + 0.5 * hdata.GetAxis(0).GetBinWidth(ptbin2)

    # print information
    print "Iteration %i, pT-bin %i, pT-range %.1f - %.1f" % (counter, ptbin1, ptmin, ptmax)

    # Create and append histogram of invariant mass
    hname = "h_invMass_%i" % (counter)
    htitle = "pT %.1f - %.1f" % (ptmin, ptmax)
    hists_inv_mass.append( hdata.Projection(1) )
    hists_inv_mass[-1].SetName(hname)
    hists_inv_mass[-1].SetTitle(htitle)

print hists_inv_mass

# create empty list to store canvases
canvases_inv_mass = []


for counter, hist in enumerate( hists_inv_mass ):

    # Create canvas and draw invariant mass plot
    cname = "iteration%i" % (counter)
    canvases_inv_mass.append( TCanvas(cname,cname) )
    hist.Draw("")

    print hist.GetTitle()


##********
#  /* set cut for direct photons */
#  h_inv_mass_allpT->GetAxis(3)->SetRange( 0 , 0 );
#
#  bini = 7;
#  for ( int ploti = 1; ploti < 17; ploti++ )
#    {
#      cout << "Plot " << ploti << " => bin " << bini << " from "
#           << h_inv_mass_allpT->GetAxis(0)->GetBinCenter(bini) - 0.5*h_inv_mass_allpT->GetAxis(0)->GetBinWidth(bini)
#           << " to "
#           << h_inv_mass_allpT->GetAxis(0)->GetBinCenter(bini) + 0.5*h_inv_mass_allpT->GetAxis(0)->GetBinWidth(bini)
#           << endl;
#
#      TString hname("invMass_pTbin");
#      hname+=bini;
#
#      h_inv_mass_allpT->GetAxis(0)->SetRange( bini , bini );
#
#      h_inv_mass_project[ploti] = (TH1F*)h_inv_mass_allpT->Projection( 1 );
#      h_inv_mass_project[ploti]->SetName(hname);
#      h_inv_mass_project[ploti]->SetTitle(hname);
#
#      bini++;
#    }
#
#  /* Duplicate plots for pT bins and do fit */
#  TH1F* h_inv_mass_project_fit[17];
#  for ( int ploti = 1; ploti < 17; ploti++ )
#    {
#      h_inv_mass_project_fit[ploti] = (TH1F*)h_inv_mass_project[ploti]->Clone();
#
#      //  h_inv_mass_project_fit[ploti]->Fit("fit_pi0","","",0.105,0.165);
#      h_inv_mass_project_fit[ploti]->Fit("fit_pi0","","",0.05,0.300);
#
#      for ( int p = 0; p < 3; p ++ )
#        fit_pi0_gauss->SetParameter(p, fit_pi0->GetParameter(p) );
#
#      cout << "Pi0 peak integral: " << fit_pi0_gauss->Integral(0.0,1.0) << endl;
#
#    }
##********


## fit pi0 peak
#fit_signal = TF1("fit_signal", "gaus(0)" )
#fit_bckgrnd  = TF1("fit_bckgrnd", "pol3(0)" )
#fit_combi = TF1("fit_combi", "gaus(0) + pol3(3)" )
#
#fit_signal.SetParameter( 0 , 4e5 )
#fit_signal.SetParameter( 1 , 0.135 )
#fit_signal.SetParameter( 2 , 0.1 )
#
#c1 = TCanvas( 'c1' , 'Fit pi0' , 200 , 10 , 700 , 500 )
#
#hdata.Projection(1).Draw()
#hdata.Projection(1).Fit('fit_signal', "", "",0.112,0.162)
#
#fit_combi.SetParameters( fit_signal.GetParameters() )
#
#hdata.Projection(1).Fit('fit_combi', "", "",0.05,0.4)
#
#
#
## fit eta peak
#fit_signal2 = TF1("fit_signal2", "gaus(0)" )
#fit_bckgrnd2  = TF1("fit_bckgrnd2", "pol3(0)" )
#fit_combi2 = TF1("fit_combi2", "gaus(0) + pol3(3)" )
#
#fit_signal2.SetParameter( 0 , 4e5 )
#fit_signal2.SetParameter( 1 , 0.55 )
#fit_signal2.SetParameter( 2 , 0.1 )
#
#c2 = TCanvas( 'c2' , 'Fit eta' , 200 , 10 , 700 , 500 )
#
#hdata.Projection(1).Draw()
#hdata.Projection(1).Fit('fit_signal2', "", "",0.5,0.6)
#
#fit_combi2.SetParameters( fit_signal2.GetParameters() )
#
#hdata.Projection(1).Fit('fit_combi2', "", "",0.4,0.7)



# close input file
f.Close()
