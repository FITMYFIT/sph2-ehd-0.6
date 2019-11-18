#deltat=0.025, output 0 2 4 6
set terminal 'wxt'
set xlabel "x"
set ylabel "{/Symbol r}_e"
set xrange [-0.51:0.51]
plot "y=0-Rhoe-BulkRelaxxxEHDBulkRelax1.dat" using 1:2 w points pt 6 ps 1.5 lw 2 lc 3 title "t=0.025 SPH", "y=0-Rhoe-BulkRelaxxxEHDBulkRelax1.dat" using 1:3 w lines lw 2 lc 7 title "t=0.025 Exact",\
     "y=0-Rhoe-BulkRelaxxxEHDBulkRelax80.dat" using 1:2 w points pt 6 ps 1.5 lw 2 lc 3 title "t=2  SPH", "y=0-Rhoe-BulkRelaxxxEHDBulkRelax80.dat" using 1:3 w lines lw 2 lc 7 title "t=2  Exact",\
     "y=0-Rhoe-BulkRelaxxxEHDBulkRelax160.dat" using 1:2 w points pt 6 ps 1.5 lw 2 lc 3 title "t=4  SPH", "y=0-Rhoe-BulkRelaxxxEHDBulkRelax160.dat" using 1:3 w lines lw 2 lc 7 title "t=4  Exact",\
     "y=0-Rhoe-BulkRelaxxxEHDBulkRelax240.dat" using 1:2 w points pt 6 ps 1.5 lw 2 lc 3 title "t=6  SPH", "y=0-Rhoe-BulkRelaxxxEHDBulkRelax240.dat" using 1:3 w lines lw 2 lc 7 title "t=6  Exact"
set terminal postscript eps color solid "Helvetica" 20 enhanced
set output "EHDBulkRelax-y=0.eps"
replot
set output
set term wxt
