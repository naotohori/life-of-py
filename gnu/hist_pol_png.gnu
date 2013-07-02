set term png enhanced
#{{no}transparent} {{no}interlace}
#{{no}truecolor} {rounded|butt}
#{tiny | small | medium | large | giant}
#{font <face> {<pointsize>}}
#{size <x>,<y>} {{no}crop}
#{{no}enhanced}
#{<color0> <color1> <color2> ...}
set size ratio 0.3
set view 180,0
#set xl "{/Symbol f}"
#set yl "{/Symbol q}" rotate by 90
set xr [180:-180]
#set yr [180:0]
set xtics 45
set ytics 30
unset ztics
#set cbtics 20
set out 'hist_pol.png'
sp 'hist_pol.gnudat' i 0 w pm3d title 'Total'
set out 'hist_pol_1.png'
sp 'hist_pol.gnudat' i 1 w pm3d title 'r < 100'
set out 'hist_pol_2.png'
sp 'hist_pol.gnudat' i 2 w pm3d title '100 < r < 150'
set out 'hist_pol_3.png'
sp 'hist_pol.gnudat' i 3 w pm3d title '150 < r'
