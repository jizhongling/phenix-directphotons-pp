import csv
import numpy as np
import matplotlib.pyplot as plt

mapfiles = [
    'warnmap-final/Warnmap_Run13pp510MinBias_Final.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_0_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_1_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_2_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_3_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange1_ERange_erange_0_to_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange2_ERange_erange_0_to_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange3_ERange_erange_0_to_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange4_ERange_erange_0_to_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange5_ERange_erange_0_to_4_nsigma10_niter10.txt'
    ]

# array of total towers per sector
ntower_total = [
    2592 ,
    2592 ,
    2592 ,
    2592 ,
    2592 ,
    2592 ,
    4608 ,
    4608
    ]

# loop over all mapfiles and load maps into arrays
warnmaps = list()

for i, mapfile in enumerate( mapfiles ):

    warnmap = np.loadtxt( mapfile, delimiter=' ' ).astype(np.int32)
    warnmaps.append( warnmap )

    print 'map: %i file: %s' % ( i, mapfiles[i] )


# store hot channel count and hot fraction for each sector / mapping file in numpy array
nmaps = len(mapfiles)
count_hot = np.zeros( (8, nmaps) )
frac_hot = np.zeros( (8, nmaps) )
frac_hot+=1

for i, warnmap in enumerate( warnmaps ):

    for sector in np.arange(0,8):

        # count hot channels
        nhot = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] == 50 ) ).sum()
        frachot = float(nhot) / float(ntower_total[sector])

        count_hot[sector][i] = nhot
        frac_hot[sector][i] = frachot

        print 'map: %i sector: %i hot fraction: %.2f  total:  %i' % ( i, sector, frachot, nhot )

    print ' '

print 'Total number of hot channels per sector (row) and warnmap file (column):'
print count_hot
print ' '

print 'Hot channel fraction per sector (row) and warnmap file (column):'
print frac_hot
print ' '

for i, mapfile in enumerate( mapfiles ):
    print 'map: %i file: %s' % ( i, mapfiles[i] )

# Plot 1: Compare different energy ranges

# first plot overall number of
plt.bar( np.arange(0,8) , count_hot[:,0], 0.95, color='crimson', label='full' )

# second plot individual energy ranges
# number of energy ranges to compare
neranges=5
widthbar=0.95/float(neranges)

plt.bar( np.arange(0,8)+0*widthbar , count_hot[:,1], widthbar, color='yellowgreen', label='0.5 GeV - 1.5 GeV' )
plt.bar( np.arange(0,8)+1*widthbar , count_hot[:,2], widthbar, color='gold', label='1.5 GeV - 3 GeV' )
plt.bar( np.arange(0,8)+2*widthbar , count_hot[:,3], widthbar, color='lightskyblue', label='3.0 GeV - 6.0 GeV' )
plt.bar( np.arange(0,8)+3*widthbar , count_hot[:,4], widthbar, color='peachpuff', label='6.0 GeV - 9.0 GeV' )
plt.bar( np.arange(0,8)+4*widthbar , count_hot[:,5], widthbar, color='mediumturquoise', label='> 9 GeV' )

plt.xlabel('sector')
plt.ylabel('# hot towers')

plt.axis([0,7.9,0,180])
plt.legend( loc='best' )

plt.savefig('plots/warnmap_stat_summary1.png')
plt.show()

# Plot 2: Compare different run ranges

# first plot overall number of
plt.bar( np.arange(0,8) , count_hot[:,0], 0.95, color='crimson', label='full' )

# second plot individual energy ranges
# number of energy ranges to compare
nrunranges=5
widthbar=0.95/float(nrunranges)

plt.bar( np.arange(0,8)+0*widthbar , count_hot[:,6], widthbar, color='yellowgreen', label='run range 1' )
plt.bar( np.arange(0,8)+1*widthbar , count_hot[:,7], widthbar, color='gold', label='run range 2' )
plt.bar( np.arange(0,8)+2*widthbar , count_hot[:,8], widthbar, color='lightskyblue', label='run range 3' )
plt.bar( np.arange(0,8)+3*widthbar , count_hot[:,9], widthbar, color='peachpuff', label='run range 4' )
plt.bar( np.arange(0,8)+4*widthbar , count_hot[:,10], widthbar, color='mediumturquoise', label='run range 5' )

plt.xlabel('sector')
plt.ylabel('# hot towers')

plt.axis([0,7.9,0,180])
plt.legend( loc='best' )

plt.savefig('plots/warnmap_stat_summary2.png')
plt.show()
