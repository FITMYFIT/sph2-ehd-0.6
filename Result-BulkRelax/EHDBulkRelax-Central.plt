set terminal 'wxt'
set xlabel "t"
set ylabel "{/Symbol r}_e"
set xrange [0:7]
set yrange [0:8]
plot "MaxRhoe-BulkRelaxxxEHDBulkRelax.dat" using 1:2 w points pt 6 ps 1.5 lc 3 lw 2 title "SPH", "MaxRhoe-BulkRelaxxxEHDBulkRelax.dat" using 1:3 w lines lw 2 lc 7  title "Exact"
set terminal postscript eps color solid linewidth 2 "Helvetica" 20 enh
set output "EHDBulkRelax-MaxRhoe.eps"
replot
set output
set term wxt