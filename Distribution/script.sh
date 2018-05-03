#!/bin/bash
clear
# NUMBER OF THE SAMPLES
N=100000

g++ Dist.cpp
./a.out ${N}

gnuplot plot.plt
