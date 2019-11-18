set terminal 'x11'
set xlabel "t"
set ylabel "Q"
set xrange [0:20]
set yrange [0:1.2]
plot "Q-evolxxEHDIsoCondCylinder.dat" using 1:2 w points pt 6 ps 1 lc 3 lw 2 title "SPH"
    
set terminal postscript eps color solid linewidth 2 "Helvetica" 20 enh
set output "EHDIsoCondCylinder-Q-evol-N128.eps"
replot
set output
set term wxt