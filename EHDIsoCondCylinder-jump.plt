set terminal 'x11'
set xlabel "r"
set ylabel "E_r"
set xrange [0:0.5]
set yrange [-0.002:0.01]
R=0.05
rhoini=0.5
S=0.009033203 #the area with R=0.05 in SPH model
epsilon2=2
E(x)=(x<R?0:(S*rhoini/(2*3.1415926*epsilon2)/x))
plot "Er-xxEHDIsoCondCylinder700.dat" using 1:2 w points pt 6 ps 1 lc 3 lw 2 title "SPH", \
     E(x) w l t "Analytical"
# "Er-xxEHDIsoCondCylinder800.dat" using 1:3 w lines  lc 3 lw 2 title "EXACT"

set terminal postscript eps color solid linewidth 2 "Helvetica" 20 enh
set output "EHDIsoCondCylinder-Q-evol.eps"
replot
set output
set term wxt