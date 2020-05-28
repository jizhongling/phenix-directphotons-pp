# Use errors from ALL results from Run9 (from AN741 and Robert Bennetts thesis) and integrated delivered luminosity values from RHIC run page
# to estimate reduction in uncertainties by using additional data sets. This assumes the uncertainties from the Run9 analysis are dominated
# by statistical uncertainties, which seems to be a valid assumption. Robert Bennetts thesis includes tables with the values for ALL and dALL.

# import packages
import matplotlib.pyplot as plt
import numpy as np
import math

import csv
from scipy.interpolate import interp1d
from scipy.interpolate import spline

# PHENIX recorded luminosities from The RHIC Spin Program: Achievements and Opportunities (1/2015) in units of pb-1
phenix_lumi_recorded = {
    'Run6pp200':7.5 ,
    'Run9pp200':16 ,
    'Run13pp510':155}

# Print used luminosity values
for i, r in enumerate(phenix_lumi_recorded):
    print "Lumi entry %d: %s -> %f pb^-1" % (i, r, phenix_lumi_recorded[r])


# ALL and dALL values from Run9 analysis according to table in Robert Bennett's thesis (correspond to plot in AN741)
l_pT_bin_min_R6 = [ 5, 6, 7,  8, 10, 12 ]
l_pT_bin_max_R6 = [ 6, 7, 8, 10, 12, 15 ]

l_ALL_R6 = [ -0.07 , 0.058 , -0.007 , -0.056 , -0.027 , -0.139 ]
l_dALL_R6 = [ 0.071 , 0.087 , 0.0996 , 0.0908 , 0.129 , 0.161 ]

l_pT_bin_ctr_R6 = [ (x + y)/2.0 for x, y in zip(l_pT_bin_min_R6, l_pT_bin_max_R6)]


# Cross section results Direct Photon production sqrt(s) = 200 GeV as function of pT bin
xsec_directphoton_pp200_pT = np.empty((0));
xsec_directphoton_pp200_xsec = np.empty((0));
with open('xsec_directphoton_pp200.txt', 'rb') as csvfile:
    linereader = csv.reader(csvfile, delimiter='\t', quotechar='/')
    for i, line in enumerate(linereader):
        if line[0].startswith('#'):
            continue
        xsec_directphoton_pp200_pT = np.append( xsec_directphoton_pp200_pT, float( line[0] ) )
        xsec_directphoton_pp200_xsec = np.append( xsec_directphoton_pp200_xsec, float( line[1] ) )
csvfile.close()

print xsec_directphoton_pp200_pT
print xsec_directphoton_pp200_xsec

# Cross section results pi0 production sqrt(s) = 200 GeV as function of pT bin
xsec_pi0_pp200_pT = np.empty((0));
xsec_pi0_pp200_xsec = np.empty((0));
with open('xsec_pi0_pp200.txt', 'rb') as csvfile:
    linereader = csv.reader(csvfile, delimiter=',', quotechar='/')
    for line in linereader:
        if line[0].startswith('#'):
            continue
        xsec_pi0_pp200_pT = np.append( xsec_pi0_pp200_pT, float( line[0] ) )
        xsec_pi0_pp200_xsec = np.append( xsec_pi0_pp200_xsec, float( line[1] ) )
csvfile.close()

print xsec_pi0_pp200_pT
print xsec_pi0_pp200_xsec

# Cross section results pi0 production sqrt(s) = 510 GeV as function of pT bin
xsec_pi0_pp510_pT = np.empty((0));
xsec_pi0_pp510_xsec = np.empty((0));
with open('xsec_pi0_pp510.txt', 'rb') as csvfile:
    linereader = csv.reader(csvfile, delimiter=',', quotechar='/')
    for line in linereader:
        if line[0].startswith('#'):
            continue
        xsec_pi0_pp510_pT = np.append( xsec_pi0_pp510_pT, float( line[0] ) )
        xsec_pi0_pp510_xsec = np.append( xsec_pi0_pp510_xsec, float( line[1] ) )
csvfile.close()

print zip( xsec_pi0_pp510_pT, xsec_pi0_pp510_xsec )

# create interpolation objects for cross section
interpol_directphoton_pp200 = interp1d(xsec_directphoton_pp200_pT, xsec_directphoton_pp200_xsec)
interpol_pi0_pp200 = interp1d(xsec_pi0_pp200_pT, xsec_pi0_pp200_xsec)
interpol_pi0_pp510 = interp1d(xsec_pi0_pp510_pT, xsec_pi0_pp510_xsec)

# estimate cross section ratio for direct photons at 510 GeV / direct photons at 200 GeV at pT bins for ALL
xsec_directphoton_pp510_over_pp200 = [ ( interpol_pi0_pp510(x) / interpol_pi0_pp200(x) )
                                       for x in l_pT_bin_ctr_R6 ]
print "Cross section ratios:"
print xsec_directphoton_pp510_over_pp200

# plot direct photon cross sections
plt.errorbar( xsec_directphoton_pp200_pT, xsec_directphoton_pp200_xsec , fmt='bs' , label='$\sqrt{s} = 200$ GeV' )
#plt.errorbar( xsec_directphoton_pp510_pT, xsec_directphoton_pp510_xsec , fmt='rs' , label='$\sqrt{s} = 510 PROJECTED$ GeV' )
plt.plot( l_pT_bin_ctr_R6 , interpol_directphoton_pp200( l_pT_bin_ctr_R6 ) , 'b-' )
plt.yscale('log')
plt.xlim([0,30])
plt.ylim([1e-2,1e4])
plt.axhline(y=0.,color='k',ls='dashed')
plt.xlabel('p$_T$ [GeV/c]')
plt.ylabel('direct photon cross section [pb*GeV^-2]$')
plt.legend(loc='best')
#plt.savefig('ALL_pT_R6.png', dpi=400, bbox_inches='tight')
plt.show()
plt.close()

# plot pi0 cross sections
plt.errorbar( xsec_pi0_pp200_pT, xsec_pi0_pp200_xsec , fmt='bs' , label='$\sqrt{s} = 200$ GeV' )
plt.errorbar( xsec_pi0_pp510_pT, xsec_pi0_pp510_xsec , fmt='rs' , label='$\sqrt{s} = 510$ GeV' )
plt.plot( l_pT_bin_ctr_R6 , interpol_pi0_pp200( l_pT_bin_ctr_R6 ) , 'b-' )
plt.plot( l_pT_bin_ctr_R6 , interpol_pi0_pp510( l_pT_bin_ctr_R6 ) , 'r-' )
plt.yscale('log')
plt.xlim([0,30])
plt.ylim([1e-10,10])
plt.axhline(y=0.,color='k',ls='dashed')
plt.xlabel('p$_T$ [GeV/c]')
plt.ylabel('pi0 cross section [mb*GeV^-2]$')
plt.legend(loc='best')
#plt.savefig('ALL_pT_R6.png', dpi=400, bbox_inches='tight')
plt.show()
plt.close()



# GRSV theory curve
grsv_pT = np.array( [ 5., 6.024, 6.183, 6.642, 6.815, 6.968, 7.107, 7.281, 7.726, 7.879, 8.039, 8.178, 8.324, 8.484, 8.617, 8.783, 8.949, 9.102, 9.284, 9.422, 9.561, 9.714, 9.867, 10.05, 10.19, 10.33, 10.49, 10.62, 10.82, 10.94, 11.1, 11.25, 11.41, 11.57, 11.72, 11.88, 12.06, 12.2, 12.34, 12.49, 12.64, 12.81, 12.96, 13.11, 13.27, 13.42, 13.57, 6.34, 6.49, 7.41, 7.58, 13.75, 13.9, 14.2, 14.37 ] )
grsv_ALL = np.array( [ 0.0031, 0.006276, 0.006693, 0.008716, 0.008287, 0.008551, 0.01001, 0.01007, 0.01113, 0.01124, 0.01147, 0.01129, 0.01253, 0.01235, 0.01274, 0.01344, 0.0137, 0.01457, 0.01558, 0.01522, 0.01523, 0.01595, 0.01611, 0.01657, 0.01676, 0.01749, 0.01746, 0.01859, 0.01866, 0.01908, 0.02039, 0.0204, 0.02111, 0.02137, 0.02138, 0.02266, 0.02258, 0.02313, 0.02377, 0.02381, 0.02503, 0.02493, 0.02526, 0.02659, 0.02663, 0.02656, 0.02709, 0.00703, 0.007994, 0.009476, 0.009963, 0.02798, 0.02798, 0.02896, 0.03088 ] )
grsv_xnew = np.linspace( grsv_pT.min(), grsv_pT.max(), 300)
grsv_smooth = spline( grsv_pT, grsv_ALL, grsv_xnew)
plt.plot( grsv_xnew, grsv_smooth, color='r', label='GRSV-std $\sqrt{s} = 200$ GeV' )

# interpolate for ALL from GRSV
interpol_l_ALL_0 = interp1d(grsv_pT, grsv_ALL)
l_ALL_0 = [ interpol_l_ALL_0(x) for x in l_pT_bin_ctr_R6 ]

# calculate projected uncertainties for other runs based on luminosity
l_dALL_projected_R13 = [ x
                         * math.sqrt( phenix_lumi_recorded[ 'Run6pp200' ] )
                         / math.sqrt( ( phenix_lumi_recorded[ 'Run13pp510' ] ) * xsec_directphoton_pp510_over_pp200[i] )
                         for i,x in enumerate(l_dALL_R6) ]

l_dALL_projected_R6_R9 = [ x * math.sqrt( phenix_lumi_recorded[ 'Run6pp200' ] ) / math.sqrt( ( phenix_lumi_recorded[ 'Run9pp200' ] +
                                                                                              phenix_lumi_recorded[ 'Run6pp200' ] ))
                           for x in l_dALL_R6 ]

l_dALL_projected_R9 = [ x * math.sqrt( phenix_lumi_recorded[ 'Run6pp200' ] ) / math.sqrt( phenix_lumi_recorded[ 'Run9pp200' ] )
                        for x in l_dALL_R6 ]

# offset x-values for better visibility
l_pT_bin_ctr_R13 = [ x + 0.3 for x in l_pT_bin_ctr_R6 ]

# offset x-values for better visibility
l_pT_bin_ctr_R9 = [ x + 0.15 for x in l_pT_bin_ctr_R6 ]
l_pT_bin_ctr_R6_R9 = [ x + 0.15 for x in l_pT_bin_ctr_R6 ]


# plot Run6 and other projected uncertainties, collection 1
plt.errorbar( l_pT_bin_ctr_R6 , l_ALL_R6 , yerr=l_dALL_R6 , fmt='sk' , capsize=3 , label='Run 6 $\sqrt{s} = 200$ GeV preliminary' )
plt.errorbar( l_pT_bin_ctr_R9 , l_ALL_0 , yerr=l_dALL_projected_R9 , fmt='sb' , capsize=3 , label='Run 9 $\sqrt{s} = 200$ GeV projected' )
plt.errorbar( l_pT_bin_ctr_R13 , l_ALL_0 , yerr=l_dALL_projected_R13 , fmt='sr' , capsize=3 , label='Run 13 $\sqrt{s} = 510$ GeV projected' )

plt.xlim([4.5,15])
plt.ylim([-0.3,0.2])
plt.axhline(y=0.,color='k',ls='dashed')
plt.xticks( np.arange(6,16,2) )
plt.yticks( np.arange(-0.3,0.2,0.1) )
plt.minorticks_on()
plt.grid()
plt.xlabel('p$_T$ [GeV/c]')
plt.ylabel('A$_{LL}$')
plt.legend(loc='lower left')
plt.savefig('projected_dALL_pT_Run13.png', dpi=400, bbox_inches='tight', transparent=False)
plt.show()
plt.close()


# DONE
