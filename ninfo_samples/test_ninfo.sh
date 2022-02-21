echo '################################### 1zih'
cd 1zih
../../ninfo_RNA.py --model DT13 --pdb 1zih.cg.pdb --hbfile 1zih.hb.list 1zih.ninfo_test
diff -sq 1zih.ninfo 1zih.ninfo_test
#rm 1zih.ninfo_test
cd ..
echo ''

echo '################################### bwyv13'
cd bwyv13
../../ninfo_RNA.py --model DT13 --pdb BWYV_coordinates_format.pdb --hbfile bwyv.hb.list bwyv_RNA13.ninfo_test
diff -sq bwyv_RNA13.ninfo bwyv_RNA13.ninfo_test
#rm bwyv_RNA13.ninfo_test
cd ..
echo ''

echo '################################### circU10'
cd circU10
../../ninfo_RNA.py --model DT13 --seq UUUUUUUUUU --circ circU10.ninfo_test
diff -sq circU10.ninfo circU10.ninfo_test
#rm circU10.ninfo_test
cd ..
echo ''

echo '################################### bwyv19'
cd bwyv19
../../ninfo_RNA.py --model NHT19 --pdb 1l2x.cg.pdb --hbfile bwyv.hb.dat --tstfile bwyv.tst.dat --end5 P --end3 P --nnhb 3 --exvfile bwyv19.exv.dat_test bwyv19.ninfo_test
diff -sq bwyv19.ninfo bwyv19.ninfo_test
#rm bwyv19.ninfo_test bwyv19.exv.dat_test
cd ..
echo ''

echo '################################### hTRHP19'
cd hTRHP19
../../ninfo_RNA.py --model NHT19 --pdb 1na2.cg.pdb --hbfile htrhp.hb.dat --end5 P --end3 P --nnhb 3 --exvfile htrhp19.exv.dat_test htrhp19.ninfo_test
diff -sq htrhp19.ninfo htrhp19.ninfo_test
#rm htrhp19.ninfo_test htrhp19.exv.dat_test
cd ..
echo ''
