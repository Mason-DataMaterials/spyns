set term png

# Add a vertical dotted line at x=0 to show centre (mean) of distribution.
set yzeroaxis

# Each bar is half the (visual) width of its x-range.
#set boxwidth 0.05 absolute
#set style fill solid 1.0 noborder

bin_width = 1.5;

bin_number(x) = floor(x/bin_width)

rounded(x) = bin_width * bin_number(x) 


set output "theta.png"
set title "theta"
plot 'output' using (rounded($2)):(2) smooth frequency with boxes notitle



