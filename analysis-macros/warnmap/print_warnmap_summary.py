import csv
import numpy as np

mapfiles = [
    'warnmap-final/Warnmap_Run13pp510_Final.txt'
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
        print 'sector: %i map: %i hot fraction: %.2f , count: %i' % ( sector, i, frachot, nhot )

    print ' '

print "*** Number of dead / uncalibrated channels per sector grouped by sector ***"
for sector in np.arange(0,8):

    for i, warnmap in enumerate( warnmaps ):

        # count bad channels
        nbad = ( ( warnmap[:,0] == sector ) & ( ( warnmap[:,3] == 100 ) | ( warnmap[:,3] == 150 ) ) ).sum()
        fracbad = float(nbad) / float(ntower_total[sector])
        print 'sector: %i map: %i dead / uncalib fraction: %.2f , count: %i' % ( sector, i, fracbad, nbad )

    print ' '

print "*** Number of 'channels around hot/dead/uncalibrated' channels per sector grouped by sector ***"
for sector in np.arange(0,8):

    for i, warnmap in enumerate( warnmaps ):

        # count excluded channels
        nbad = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] == 40 ) ).sum()
        fracbad = float(nbad) / float(ntower_total[sector])
        print 'sector: %i map: %i \'excluded around bad\' fraction: %.2f , count: %i' % ( sector, i, fracbad, nbad )

    print ' '

print "*** Number of 'fiducial cut sector edge' channels per sector grouped by sector ***"
for sector in np.arange(0,8):

    for i, warnmap in enumerate( warnmaps ):

        # count excluded channels
        nbad = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] == 20 ) ).sum()
        fracbad = float(nbad) / float(ntower_total[sector])
        print 'sector: %i map: %i \'excluded at sector edge\' fraction: %.2f , count: %i' % ( sector, i, fracbad, nbad )

    print ' '

print "*** Number of live channels per sector grouped by sector ***"
for sector in np.arange(0,8):

    for i, warnmap in enumerate( warnmaps ):

        # count live channels
        nlive = ( ( warnmap[:,0] == sector ) & ( warnmap[:,3] < 20 ) ).sum()
        fraclive = float(nlive) / float(ntower_total[sector])
        print 'sector: %i map: %i live fraction: %.3f total: %i' % ( sector, i, fraclive, nlive )

    print ' '
