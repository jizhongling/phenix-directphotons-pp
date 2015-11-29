import csv
import numpy as np

mapfiles = [
    'warnmap-output/Warnmap_Run13pp510MinBias_AND_erange_0_to_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_0_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_1_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_2_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_3_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange1_MergeEnergyRanges_erange_0_to_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange2_MergeEnergyRanges_erange_0_to_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange3_MergeEnergyRanges_erange_0_to_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange4_MergeEnergyRanges_erange_0_to_4_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_RunRange5_MergeEnergyRanges_erange_0_to_4_nsigma10_niter10.txt'
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

print "*** Number of hot channels per sector grouped by warnmap file ***"
for i, warnmap in enumerate( warnmaps ):

    for sector in np.arange(0,8):
        # count hot channels
        nhot = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] == 50 ) ).sum()
        frachot = float(nhot) / float(ntower_total[sector])
        print 'map: %i sector: %i hot fraction: %.2f  total:  %i' % ( i, sector, frachot, nhot )

    print ' '

print "*** Number of hot channels per sector grouped by sector ***"
for sector in np.arange(0,8):

    for i, warnmap in enumerate( warnmaps ):

        # count hot channels
        nhot = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] == 50 ) ).sum()
        frachot = float(nhot) / float(ntower_total[sector])
        print 'sector: %i map: %i hot fraction: %.2f total: %i' % ( sector, i, frachot, nhot )

    print ' '

print "*** Number of live channels per sector grouped by sector ***"
for sector in np.arange(0,8):

    for i, warnmap in enumerate( warnmaps ):

        # count live channels
        nlive = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] < 20 ) ).sum()
        fraclive = float(nlive) / float(ntower_total[sector])
        print 'sector: %i map: %i live fraction: %.3f total: %i' % ( sector, i, fraclive, nlive )

    print ' '
