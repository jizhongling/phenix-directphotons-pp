# Use errors from ALL results from Run9 (from AN741 and Robert Bennetts thesis) and integrated delivered luminosity values from RHIC run page
# to estimate reduction in uncertainties by using additional data sets. This assumes the uncertainties from the Run9 analysis are dominated
# by statistical uncertainties, which seems to be a valid assumption. Robert Bennetts thesis includes tables with the values for ALL and dALL.

# import packages
import matplotlib.pyplot as plt
import numpy as np
import math

# RHIC delivered luminosity from RHIC run page http://www.rhichome.bnl.gov/RHIC/Runs/ in units of pb-1
rhic_lumi_delivered = {
    'Run6pp200':88.6 ,
    'Run9pp200':114 ,
    'Run12pp200':74.0 ,
    'Run15pp200':382 }

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


# set ALL values to 0 to focus on uncertainties
l_ALL_0 = [ 0.0 * x for x in l_ALL_R6 ]


# calculate projected uncertainties for other runs based on luminosity
l_dALL_projected_R6_R9 = [ x * math.sqrt( phenix_lumi_recorded[ 'Run6pp200' ] ) / math.sqrt( ( phenix_lumi_recorded[ 'Run9pp200' ] +
                                                                                              phenix_lumi_recorded[ 'Run6pp200' ] ))
                           for x in l_dALL_R6 ]

l_dALL_projected_R9 = [ x * math.sqrt( phenix_lumi_recorded[ 'Run6pp200' ] ) / math.sqrt( phenix_lumi_recorded[ 'Run9pp200' ] )
                        for x in l_dALL_R6 ]

# offset x-values for better visibility
l_pT_bin_ctr_R9 = [ x + 0.15 for x in l_pT_bin_ctr_R6 ]
l_pT_bin_ctr_R6_R9 = [ x + 0.15 for x in l_pT_bin_ctr_R6 ]


# Plot ALL vs pT from Run 6 (Rober Bennet thesis and AN741
fig_ALL_pT = plt.figure()
plt.errorbar( l_pT_bin_ctr_R6 ,l_ALL_R6 , yerr=l_dALL_R6 , fmt='s' , label='Run 6 $\sqrt{s} = 200$ GeV' )
plt.xlim([4.5,15])
plt.ylim([-0.5,0.5])
plt.axhline(y=0.,color='k',ls='dashed')
plt.xticks( np.arange(6,16,2) )
plt.yticks( np.arange(-0.4,0.6,0.2) )
plt.minorticks_on()
plt.grid()
plt.xlabel('p$_T$ [GeV/c]')
plt.ylabel('A$_{LL}$')
plt.legend(loc='best')
plt.savefig('ALL_pT_R6.png', dpi=400, bbox_inches='tight')
#plt.show()
plt.close()


# plot Run6 and other projected uncertainties, collection 1
fig_projected_dALL_pT = plt.figure()
plt.errorbar( l_pT_bin_ctr_R6 , l_ALL_0 , yerr=l_dALL_R6 , fmt='sk' , label='Run 6 $\sqrt{s} = 200$ GeV' )
plt.errorbar( l_pT_bin_ctr_R9 , l_ALL_0 , yerr=l_dALL_projected_R9 , fmt='sb' , label='Run 9 $\sqrt{s} = 200$ GeV' )

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
plt.savefig('projected_dALL_pT.png', dpi=400, bbox_inches='tight', transparent=True)
plt.show()
plt.close()



# plot Run6 and other projected uncertainties, collection 2
fig_projected_dALL_pT = plt.figure()
plt.errorbar( l_pT_bin_ctr_R6 , l_ALL_0 , yerr=l_dALL_R6 , fmt='sk' , label='Run 6 $\sqrt{s} = 200$ GeV' )
plt.errorbar( l_pT_bin_ctr_R6_R9 , l_ALL_0 , yerr=l_dALL_projected_R6_R9 , fmt='sb' , label='Run 6+9 $\sqrt{s} = 200$ GeV' )

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
plt.savefig('projected_dALL_pT_v2.png', dpi=400, bbox_inches='tight')
plt.show()
#plt.close()

# DONE
