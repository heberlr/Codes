clear
reset
set encoding utf8
set term postscript eps enhanced color
set key off
set border 3
set output 'Gaussian.eps'
set xrange [0:1]

# Add a vertical dotted line at x=0 to show centre (mean) of distribution.
set yzeroaxis

# Each bar is half the (visual) width of its x-range.
set boxwidth 0.01 absolute
set style fill solid 1.0 noborder

bin_width = 0.01;

bin_number(x) = floor(x/bin_width)

rounded(x) = bin_width * ( bin_number(x) + 0.5 )

plot 'gaussian.dat' using (rounded($1)):(1) smooth frequency with boxes
set output 'Uniform.eps'
plot 'uniform.dat' using (rounded($1)):(1) smooth frequency with boxes
