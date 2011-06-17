#!/bin/bash

for i in `ls`
do
    test -d $i && make -C $i clean
done