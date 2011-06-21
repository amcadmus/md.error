#!/bin/bash

l1=70
l2=80
g1=20
g2=130

cd tools/
svn up
cd ..
make -C tools/cal.density.profile -j4

tools/cal.density.profile/density --check-point-liquid-1 $l1 --check-point-liquid-2 $l2 --check-point-gas-1 $g1 --check-point-gas-2 $g2 -s 5000