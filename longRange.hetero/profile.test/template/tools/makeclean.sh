#!/bin/bash

for i in `ls`
do
    test -d $i && test -f $i/Makefile && make -C $i clean
done
