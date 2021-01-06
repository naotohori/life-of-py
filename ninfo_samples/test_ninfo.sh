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

cd bwyv19
../../ninfo_RNA.py --model NHT19 --pdb 1l2x.cg.pdb --hbfile bwyv.hb.dat --tstfile bwyv.tst.dat --end5 P --end3 P --nnhb 3 --exvfile bwyv19.exv.dat_test bwyv19.ninfo_test
diff -q bwyv19.ninfo bwyv19.ninfo_test
rm bwyv19.ninfo_test bwyv19.exv.dat_test
cd ..
