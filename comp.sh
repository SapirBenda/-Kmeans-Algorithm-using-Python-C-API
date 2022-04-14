#!/bin/bash
# Script to compile and execute a c program
gcc -ansi -Wall -Wextra -Werror -pedantic-errors common.c kmeans.c spkmeans.c -lm -o spkmeans


./spkmeans 0 wam test.csv
./spkmeans 0 ddg test.csv
./spkmeans 0 lnorm test.csv
./spkmeans 0 jacobi Lnorm.txt
./spkmeans 0 jacobi test.csv
