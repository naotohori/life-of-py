set term postscript enhanced color 20 lw 3 solid
set out 'hist_pol.ps'
set size ratio 0.3
set view 180,0
set xl "{/Symbol f}"
set yl "{/Symbol q}" rotate by 90
set xr [180:-180]
#set yr [180:0]
set xtics 45
set ytics 30
unset ztics
#set cbtics 20
sp 'hist_pol.gnudat' i 0 w pm3d title 'Total'
sp 'hist_pol.gnudat' i 1 w pm3d title 'r < 100'
sp 'hist_pol.gnudat' i 2 w pm3d title '100 < r < 150'
sp 'hist_pol.gnudat' i 3 w pm3d title '150 < r'
