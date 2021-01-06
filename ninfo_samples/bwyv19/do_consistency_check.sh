
# At semiexplicit/examples/
../semiexplicit -p 1l2x.pdb -b hbond.dat -k stack.dat -m $bin/maxi_explicit -u $bin/uvv/pmf_Mg_P_1264 -T 80 -M 0.000 -K 0.5 -s 10000 -v "700 700 700" -a 1 -e 1 -C 30 -S 111 -z 29 -d 10 1> run.out 2> run.err

cp hbond.log ~/python/cafysis/ninfo_samples/bwyv19/semiexplicit.hbond.log

./ninfo_hbond_check_consistency.py semiexplicit.hbond.log bwyv19.ninfo
