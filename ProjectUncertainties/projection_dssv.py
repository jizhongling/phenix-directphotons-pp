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



# Change pT to xT
l_xT_bin_ctr_R6 = [ x / 100. for x in l_pT_bin_ctr_R6 ]
l_xT_bin_ctr_R9 = [ x / 100. for x in l_pT_bin_ctr_R6 ]
l_xT_bin_ctr_R13 = [ x / 255. for x in l_pT_bin_ctr_R6 ]

# DSSV theory curve 200 GeV
dssv200_pT = np.array( [ 3.000000E+00, 4.000000E+00, 5.000000E+00, 6.000000E+00, 7.000000E+00, 8.000000E+00, 9.000000E+00, 1.000000E+01, 1.100000E+01, 1.200000E+01, 1.300000E+01, 1.400000E+01, 1.500000E+01, 1.600000E+01, 1.700000E+01, 1.800000E+01, 1.900000E+01, 2.000000E+01, 2.100000E+01, 2.200000E+01, 2.300000E+01, 2.400000E+01, 2.500000E+01, 2.600000E+01, 2.700000E+01, 2.800000E+01, 2.900000E+01, 3.000000E+01, 3.100000E+01, 3.200000E+01, 3.400000E+01, 3.600000E+01, 3.800000E+01, 4.000000E+01, 4.200000E+01, 4.400000E+01, 4.600000E+01, 4.800000E+01, 5.000000E+01, 5.200000E+01, 5.400000E+01, 5.600000E+01, 5.800000E+01, 6.000000E+01, 6.200000E+01, 6.400000E+01 ] )
dssv200_ALL = np.array( [ 4.439147E-04, 1.171262E-03, 2.220237E-03, 3.601331E-03, 5.306280E-03, 7.320380E-03, 9.623244E-03, 1.219952E-02, 1.502220E-02, 1.805147E-02, 2.131348E-02, 2.472844E-02, 2.828813E-02, 3.196382E-02, 3.571679E-02, 3.957304E-02, 4.352939E-02, 4.748573E-02, 5.148410E-02, 5.537836E-02, 5.929806E-02, 6.310377E-02, 6.690826E-02, 7.058333E-02, 7.413055E-02, 7.747778E-02, 8.068857E-02, 8.364492E-02, 8.645996E-02, 8.912327E-02, 9.382894E-02, 9.753721E-02, 1.002440E-01, 1.017501E-01, 1.022255E-01, 1.018049E-01, 1.006426E-01, 9.859850E-02, 9.592489E-02, 9.260725E-02, 8.887382E-02, 8.474397E-02, 8.063128E-02, 7.638462E-02, 7.222191E-02, 6.815836E-02 ] )
dssv200_xT = [ x / 100. for x in dssv200_pT ]
dssv200_xnew = np.linspace( 5./100., 15./100., 300)
dssv200_smooth = spline( dssv200_xT, dssv200_ALL, dssv200_xnew)
plt.plot( dssv200_xnew, dssv200_smooth, color='k', label='DSSV14 $\sqrt{s} = 200$ GeV' )

# interpolate for ALL from DSSV
interpol_l_ALL_200 = interp1d(dssv200_xT, dssv200_ALL)
l_ALL_200 = [ interpol_l_ALL_200(x) for x in l_xT_bin_ctr_R6 ]

# DSSV theory curve 500 GeV
dssv500_pT = np.array( [ 3.000000E+00, 4.000000E+00, 5.000000E+00, 6.000000E+00, 7.000000E+00, 8.000000E+00, 9.000000E+00, 1.000000E+01, 1.100000E+01, 1.200000E+01, 1.300000E+01, 1.400000E+01, 1.500000E+01, 1.600000E+01, 1.700000E+01, 1.800000E+01, 1.900000E+01, 2.000000E+01, 2.100000E+01, 2.200000E+01, 2.300000E+01, 2.400000E+01, 2.500000E+01, 2.600000E+01, 2.700000E+01, 2.800000E+01, 2.900000E+01, 3.000000E+01, 3.100000E+01, 3.200000E+01, 3.500000E+01, 4.000000E+01, 4.500000E+01, 5.000000E+01, 5.500000E+01, 6.000000E+01, 6.500000E+01, 7.000000E+01, 7.500000E+01, 8.000000E+01, 8.500000E+01, 9.000000E+01, 9.500000E+01, 1.000000E+02, 1.050000E+02, 1.100000E+02, 1.150000E+02, 1.200000E+02, 1.250000E+02, 1.300000E+02, 1.350000E+02, 1.400000E+02, 1.450000E+02, 1.500000E+02, 1.550000E+02, 1.600000E+02 ] )
dssv500_ALL = np.array( [ 7.387409E-05, 1.888835E-04, 3.510941E-04, 5.658036E-04, 8.355659E-04, 1.162648E-03, 1.543574E-03, 1.984241E-03, 2.477886E-03, 3.026275E-03, 3.637747E-03, 4.298877E-03, 5.014310E-03, 5.781770E-03, 6.598182E-03, 7.466651E-03, 8.382832E-03, 9.336199E-03, 1.035077E-02, 1.138738E-02, 1.248215E-02, 1.361322E-02, 1.478343E-02, 1.597633E-02, 1.720694E-02, 1.846706E-02, 1.976279E-02, 2.107464E-02, 2.242041E-02, 2.381077E-02, 2.809116E-02, 3.554049E-02, 4.329603E-02, 5.113261E-02, 5.883241E-02, 6.633830E-02, 7.354276E-02, 8.012240E-02, 8.596952E-02, 9.115242E-02, 9.542684E-02, 9.877761E-02, 1.012237E-01, 1.025992E-01, 1.033900E-01, 1.031852E-01, 1.022967E-01, 1.006371E-01, 9.844477E-02, 9.568568E-02, 9.259108E-02, 8.931748E-02, 8.589496E-02, 8.225350E-02, 7.867339E-02, 7.514518E-02 ] )
dssv500_xT = [ x / 250. for x in dssv500_pT ]
dssv500_xnew = np.linspace( 5./250., 15./250., 300)
dssv500_smooth = spline( dssv500_xT, dssv500_ALL, dssv500_xnew)
plt.plot( dssv500_xnew, dssv500_smooth, color='r', label='DSSV14 $\sqrt{s} = 500$ GeV' )

# interpolate for ALL from DSSV
interpol_l_ALL_500 = interp1d(dssv500_xT, dssv500_ALL)
l_ALL_500 = [ interpol_l_ALL_500(x) for x in l_xT_bin_ctr_R13 ]

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



# plot Run6 and other projected uncertainties, collection 1
plt.errorbar( l_xT_bin_ctr_R6 , l_ALL_R6 , yerr=l_dALL_R6 , fmt='sk' , capsize=3 , label='Run 6 $\sqrt{s} = 200$ GeV preliminary' )
#plt.errorbar( l_xT_bin_ctr_R9 , l_ALL_200 , yerr=l_dALL_projected_R9 , fmt='sb' , capsize=3 , label='Run 9 $\sqrt{s} = 200$ GeV projected' )
plt.errorbar( l_xT_bin_ctr_R13 , l_ALL_500 , yerr=l_dALL_projected_R13 , fmt='sr' , capsize=3 , label='Run 13 $\sqrt{s} = 510$ GeV projected' )

plt.xlim([0.,0.15])
plt.ylim([-0.3,0.2])
plt.axhline(y=0.,color='k',ls='dashed')
plt.xticks( np.arange(0.,0.15,0.05) )
plt.yticks( np.arange(-0.3,0.2,0.1) )
plt.minorticks_on()
plt.grid()
plt.xlabel('x$_T$')
plt.ylabel('A$_{LL}$')
plt.legend(loc='lower left')
plt.savefig('projected_dALL_pT_Run13.png', dpi=400, bbox_inches='tight', transparent=False)
plt.show()
plt.close()


# DONE
