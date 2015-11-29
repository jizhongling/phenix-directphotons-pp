import csv
import numpy as np

mapfiles1 = [
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_0_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_1_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_2_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_3_nsigma10_niter10.txt' ,
    'warnmap-output/Warnmap_Run13pp510MinBias_erange_4_nsigma10_niter10.txt'
    ]

# select set of files to merge
mapfiles=mapfiles1

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


# merge warnmap files
print "*** Number of hot channels per sector for merged warnmap ***"

warnmap_merged = warnmaps[0].copy()

for warnmap in warnmaps:
    warnmap_merged = np.maximum( warnmap_merged, warnmap )

# count hot channels
for sector in np.arange(0,8):
    nhot = ( ( warnmap_merged[:,0] == sector ) & ( warnmap_merged[:,3] == 50 ) ).sum()
    frachot = float(nhot) / float(ntower_total[sector])
    print 'sector: %i hot fraction: %.2f  total:  %i' % ( sector, frachot, nhot )

# save warnmap to txt file
np.savetxt( 'warnmap_merged_python.txt', warnmap_merged, delimiter=' ', fmt='%1.1d' )

