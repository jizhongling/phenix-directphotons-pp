import csv
import numpy as np
import matplotlib.pyplot as plt
import sys

# change font used by matplotlib
#font = { 'family' : 'normal',
#         'weight' : 'none',
#         'size'   : 18}
#
#plt.rc('font', **font)

#plt.rcParams.update({'axes.titlesize' : 16})
plt.rcParams.update({'axes.labelsize' : 18})
plt.rcParams.update({'xtick.labelsize' : 18})
plt.rcParams.update({'ytick.labelsize' : 18})

#mapfiles = ['warnmap-final/Warnmap_Run13pp510MinBias_Final.txt' ,
#            'warnmap-output/Warnmap_Run13pp510MinBias_erange_0_nsigma10_niter10.txt'
#            ]
#plotname = 'plots/warnmap_compare_Run13pp510MinBias_final_erange0.png'

#mapfiles = ['warnmap-final/Warnmap_Run13pp510MinBias_Final.txt' ,
#            'warnmap-final/Warnmap_Run13pp510ERT_Final.txt'
#            ]
#plotname = 'plots/warnmap_compare_Run13pp510MinBias_Final_Run13pp510ERT_Final.png'

#mapfiles = [ 'warnmap-output/Warnmap_Run13pp510MinBias_erange_1_nsigma10_niter10.txt' ,
#             'warnmap-output/Warnmap_Run13pp510MinBias_erange_2_nsigma10_niter10.txt'
#            ]
#plotname = 'plots/warnmap_compare_Run13pp510MinBias_erange_1_vs_2.png'

#mapfiles = [ 'warnmap-output/Warnmap_Run13pp510MinBias_erange_0_nsigma10_niter10.txt' ,
#             'warnmap-output/Warnmap_Run13pp510MinBias_erange_4_nsigma10_niter10.txt'
#            ]
#plotname = 'plots/warnmap_compare_Run13pp510MinBias_erange_0_vs_4.png'

#mapfiles = ['warnmap-final/Warnmap_Run13pp510MinBias_Final.txt' ,
#            'warnmap-final/Warnmap_Run9pp500MinBias_Final.txt'
#            ]
#plotname = 'plots/warnmap_compare_Run13pp510MinBias_Run9pp500MinBias.png'

mapfiles = [ 'warnmap-final/Warnmap_Run9pp500MinBias_Final.txt',
             'warnmap-paul/iter10_rms10/warn_Run9pp500MinBias_newFormat.txt'
            ]
plotname = 'plots/warnmap_compare_Run9pp500MinBias_me_vs_paul.png'

#mapfiles = ['warnmap-final/Warnmap_Run13pp510MinBias_Final.txt' ,
#            'warnmap-paul/iter10_rms10/warn_Run9pp500MinBias.txt'
#            ]
#plotname = 'plots/warnmap_compare_Run13pp510MinBias_test.png'


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

print 'Comparing map 0 (%s)' % ( mapfiles[0] )
print '       to map 1 (%s)' % ( mapfiles[1] )
print ' '

for i, mapfile in enumerate( mapfiles ):

    warnmap = np.loadtxt( mapfile, delimiter=' ' ).astype(np.int32)
    warnmaps.append( warnmap )

# store hot channel count and hot fraction for each sector / mapping file in numpy array
nmaps = len(mapfiles)
count_hot = np.zeros( (8, nmaps) )
count_hot_0_AND_1 = np.zeros( (8) )
count_hot_0_AND_NOT_1 = np.zeros( (8) )
count_hot_NOT_0_AND_1 = np.zeros( (8) )

#frac_hot = np.zeros( (8, nmaps) )
#frac_hot+=1

for sector in np.arange(0,8):

    for i, warnmap in enumerate( warnmaps ):
        # count hot channels
        nhot = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] == 50 ) ).sum()
        #frachot = float(nhot) / float(ntower_total[sector])

        count_hot[sector][i] = nhot
        #frac_hot[sector][i] = frachot

    nhot_0_and_1 = ( ( warnmaps[0][:,0] == sector ) &
                     ( warnmaps[0][:,0] == warnmaps[1][:,0] ) &
                     ( warnmaps[0][:,1] == warnmaps[1][:,1] ) &
                     ( warnmaps[0][:,2] == warnmaps[1][:,2] ) &
                     ( warnmaps[0][:,3] == 50 ) &
                     ( warnmaps[1][:,3] == 50 ) ).sum()

    nhot_0_and_not_1 = ( ( warnmaps[0][:,0] == sector ) &
                         ( warnmaps[0][:,0] == warnmaps[1][:,0] ) &
                         ( warnmaps[0][:,1] == warnmaps[1][:,1] ) &
                         ( warnmaps[0][:,2] == warnmaps[1][:,2] ) &
                         ( warnmaps[0][:,3] == 50 ) &
                         ( warnmaps[1][:,3] != 50 ) ).sum()

    nhot_not_0_and_1 = ( ( warnmaps[0][:,0] == sector ) &
                         ( warnmaps[0][:,0] == warnmaps[1][:,0] ) &
                         ( warnmaps[0][:,1] == warnmaps[1][:,1] ) &
                         ( warnmaps[0][:,2] == warnmaps[1][:,2] ) &
                         ( warnmaps[0][:,3] != 50 ) &
                         ( warnmaps[1][:,3] == 50 ) ).sum()

    count_hot_0_AND_1[sector] = nhot_0_and_1
    count_hot_0_AND_NOT_1[sector] = nhot_0_and_not_1
    count_hot_NOT_0_AND_1[sector] = nhot_not_0_and_1

    print 'Hot in map0                 (sector %i): %i' % (sector, count_hot[sector,0])
    print 'Hot in map1                 (sector %i): %i' % (sector, count_hot[sector,1])
    print 'Hot in map0 AND in map1     (sector %i): %i' % (sector, nhot_0_and_1)
    print 'Hot in map0 AND NOT in map1 (sector %i): %i' % (sector, nhot_0_and_not_1)
    print 'Hot NOT in map0 AND in map1 (sector %i): %i' % (sector, nhot_not_0_and_1)
    print ' '



# Plot 1: Compare overlap in warnmaps
widthbar=0.45

plt.bar( np.arange(0,8)+0*widthbar , count_hot[:,0], widthbar, color='crimson', label='Hot in map A')
plt.bar( np.arange(0,8)+1*widthbar , count_hot[:,1], widthbar, color='yellowgreen', label='Hot in map B' )
plt.bar( np.arange(0,8) , count_hot_0_AND_1[:], 2*widthbar, color='gold', label='Hot in map A && map B' , edgecolor='black' , hatch='xx' , fill=False )

plt.xlabel('sector')
plt.ylabel('# hot towers')

plt.axis([0,7.9,0,180])
plt.legend( loc='best' )

plt.savefig(plotname)
plt.show()
