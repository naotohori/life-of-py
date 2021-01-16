pdb_aa2cg.py 1na2.pdb Y 1na2.cg.pdb

ninfo_RNA.py --model NHT19 --pdb 1na2.cg.pdb --hbfile htrhp.hb.dat --end5 P --end3 P --nnhb 3 --exvfile htrhp19.exv.dat htrhp19.ninfo

../../ninfo_hbond_check_consistency.py semiexplicit.hbond.log htrhp19.ninfo
