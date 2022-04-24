#!/bin/bash
# Script to compile and execute a c program
gcc -ansi -Wall -Wextra -Werror -pedantic-errors common.c kmeans.c spkmeans.c -lm -o spkmeans

# run commands for example
./spkmeans 0 wam inputs_examples/data_points/input_1.txt
./spkmeans 0 ddg inputs_examples/data_points/input_1.txt
./spkmeans 0 lnorm inputs_examples/data_points/input_1.txt
./spkmeans 0 jacobi inputs_examples/sym_matrix/sym_matrix_input_1.txt
