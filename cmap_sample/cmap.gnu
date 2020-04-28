set view map
set size ratio 1

set zr [0:3]
set cbrange [0:3]
set palette defined (0 "white", 1 "black", 2 "blue", 3"red")


unset key
set xtics 10 out
set ytics 10 out
unset ztics
set grid


## protein L
set xr [0:73]
set yr [0:73]

set term postscript color 20 lw 2 solid 
set out 'proteinL_cmap.ps'
sp 'proteinL.gnudat' u ($1+0.5):($2+0.5):3  w pm3d
#,  'proteinL.gnudat' u ($1+0.5):($2+0.5):3 w l lw 0.1 lt 0

set term svg
set out 'proteinL_cmap.svg'
replot


