import csv
import numpy as np

mapfiles = [
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange0.txt' ,
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange1.txt' ,
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange2.txt' ,
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange3.txt' ,
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange4.txt'
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

for mapfile in mapfiles:

    warnmap = np.loadtxt( mapfile, delimiter=' ' ).astype(np.int32)
    warnmaps.append( warnmap )

print "*** Number of hot channels per sector grouped by warnmap file ***"
for i, warnmap in enumerate( warnmaps ):

    for sector in np.arange(0,8):
        # count hot channels
        nhot = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] == 50 ) ).sum()
        frachot = float(nhot) / float(ntower_total[sector])
        print 'map: %s sector: %i hot fraction: %.2f  total:  %i' % ( mapfiles[i][-11:-4], sector, frachot, nhot )

    print ' '

print "*** Number of hot channels per sector grouped by sector ***"
for sector in np.arange(0,8):

    for i, warnmap in enumerate( warnmaps ):

        # count hot channels
        nhot = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] == 50 ) ).sum()
        frachot = float(nhot) / float(ntower_total[sector])
        print 'sector: %i map: %s hot fraction: %.2f total: %i' % ( sector, mapfiles[i][-11:-4], frachot, nhot )

    print ' '
