#!/bin/bash
clear
rm *~ *# &> /dev/null
r=${1:-'100'}
p=${2:-'4'}
Morris=${3:-'0'}
K=${4:-'6'}
g++ -O3 codigo.cpp mgmres.cpp 
time ./a.out ${r} ${p} ${Morris} 7
rm *.out &> /dev/null
gnuplot saida.plt
