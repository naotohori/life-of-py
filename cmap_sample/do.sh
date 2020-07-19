./convert.py proteinL.dat 1k53_A.dssp_struct proteinL.gnudat

gnuplot cmap.gnu

sed -i -e "/<polygon fill = 'rgb(255, 255, 255)' points =/d" proteinL_cmap.svg

rm *.svg-e
