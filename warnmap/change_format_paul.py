import numpy as np
import sys

# format Paul's warnmap:
# towerid status

# output format:
# sector y-bin z-bin status

# define function to decode Paul's tower index into sector, y-bin, z-bin
def decode_index( towerid , status):#, sector , y , z ):
    itwr=0

    sector = -1
    z = -1
    y = -1

    if( towerid < 0 or towerid > 24767 ):
        print 'BAD tower id ', towerid

    else:

        if( towerid < 15552 ):
            # pbsc
            sector = towerid // 2592
            itwr = towerid % 2592
            z = itwr % 72
            y = itwr // 72
        else:
            # pbgl
            sector = 6 + ( towerid - 15552 ) // 4608
            itwr = ( towerid - 15552 ) % 4608
            z = itwr % 96
            y = itwr // 96

    #print 'Index %f -> Sector %f, y-bin %f, z-bin %f' % ( towerid, sector, y, z )
    return (sector,y,z,status)

print 'Running macro: ', sys.argv[0]
print 'This is the number of arguments: ', len(sys.argv)

if ( len(sys.argv) < 3 ):
    print 'Usage: ', sys.argv[0], '<map_input.txt> <map_output.txt>'
    sys.exit()

map_input = np.loadtxt(sys.argv[1])

print map_input
print 'Size of array from input file: ', map_input.shape

# convert index column to position column with sector, y-bin, z-bin; keep status column
decode_index_vec = np.vectorize(decode_index)

map_output = np.transpose( np.asarray( decode_index_vec( map_input[:,0],map_input[:,1]  ) ) )

# convert status definition
# status '50' to '40' (Paul 50 = around hot tower and on sector edges)
map_output[ map_output[:,3] == 50, 3 ] = 40

# status '100' to '50' (Paul 100 = hot)
map_output[ map_output[:,3] == 100, 3 ] = 50

print map_output
print 'Size of array after format changes: ', map_output.shape

np.savetxt( sys.argv[2], map_output, fmt='%i' )


