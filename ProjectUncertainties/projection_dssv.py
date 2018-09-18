# Use errors from ALL results from Run9 (from AN741 and Robert Bennetts thesis) and integrated delivered luminosity values from RHIC run page
# to estimate reduction in uncertainties by using additional data sets. This assumes the uncertainties from the Run9 analysis are dominated
# by statistical uncertainties, which seems to be a valid assumption. Robert Bennetts thesis includes tables with the values for ALL and dALL.

# import packages
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.interpolate import spline

import csv
from scipy.interpolate import interp1d

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

l_ALL_R6 = [ -0.070 , 0.058 , -0.006 , -0.056 , -0.027 , -0.139 ]
l_dALL_R6 = [ 0.071 , 0.087 , 0.100 , 0.091 , 0.113 , 0.161 ]

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



# set ALL values to 0 to focus on uncertainties
#l_ALL_0 = [ 0.0 * x for x in l_ALL_R6 ]
l_ALL_0 = [ 0.00348347, 0.00727476, 0.0101439, 0.0139859, 0.0197229, 0.0282783 ]

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
fig_projected_dALL_pT = plt.figure()
plt.errorbar( l_pT_bin_ctr_R6 , l_ALL_0 , yerr=l_dALL_R6 , fmt='sk' , label='Run 6 $\sqrt{s} = 200$ GeV' )
plt.errorbar( l_pT_bin_ctr_R9 , l_ALL_0 , yerr=l_dALL_projected_R9 , fmt='sb' , label='Run 9 $\sqrt{s} = 200$ GeV' )
plt.errorbar( l_pT_bin_ctr_R13 , l_ALL_0 , yerr=l_dALL_projected_R13 , fmt='sr' , label='Run 13 $\sqrt{s} = 510$ GeV' )

# DSSV theory curve
dssv_pT = np.array( [ 5.50388, 6.03711, 6.50681, 6.86229, 7.50971, 8.27140, 9.00771, 9.65518, 10.4930, 11.0008, 11.7244, 12.5877, 13.5018, 14.3523, 14.7459 ] ) 
dssv_ALL = np.array( [ 0.00348347, 0.00630450, 0.00727476, 0.00915544, 0.0101439, 0.0120662, 0.0139859, 0.0168186, 0.0187488, 0.0197229, 0.0225634, 0.0244962, 0.0282783, 0.0302097, 0.0320943 ] )
dssv_xnew = np.linspace( dssv_pT.min(), dssv_pT.max(), 300)
dssv_smooth = spline( dssv_pT, dssv_ALL, dssv_xnew)
plt.plot( dssv_xnew, dssv_smooth, color='r', label='GRSV-std for direct photon' )

plt.xlim([4.5,15])
plt.ylim([-0.2,0.2])
plt.axhline(y=0.,color='k',ls='dashed')
plt.xticks( np.arange(6,16,2) )
plt.yticks( np.arange(-0.2,0.2,0.05) )
plt.minorticks_on()
plt.grid()
plt.xlabel('p$_T$ [GeV/c]')
plt.ylabel('A$_{LL}$')
plt.legend(loc='best')
plt.savefig('projected_dALL_pT_Run13.png', dpi=400, bbox_inches='tight', transparent=False)
plt.show()
plt.close()


# DONE
