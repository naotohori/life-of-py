cd bwyv13
../../ninfo_RNA.py --model DT13 --pdb BWYV_coordinates_format.pdb --hbfile bwyv.hb.list bwyv_RNA13.ninfo_test
diff -q bwyv_RNA13.ninfo bwyv_RNA13.ninfo_test
rm bwyv_RNA13.ninfo_test
cd ..

cd circU10
../../ninfo_RNA.py --model DT13 --seq UUUUUUUUUU --circ circU10.ninfo_test
diff -q circU10.ninfo circU10.ninfo_test
rm circU10.ninfo_test
cd ..
