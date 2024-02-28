#!/bin/bash
# converts MNI coordinates in mm to the voxel coordinates
 
x=$(echo "($1 * -1 + 90)/2" | bc)
y=$(echo "($2 * 1 + 126)/2" | bc)
z=$(echo "($3 * 1 + 72)/2" | bc)

echo $x $y $z