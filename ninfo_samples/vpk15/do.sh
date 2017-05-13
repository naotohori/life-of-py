pdb_detect_basestack.py VPK_0_34.pdb > pdb_detect_basestack.out

cp pdb_detect_basestack.out VPK.stack.list
echo 'F  -15  -29' >> VPK.stack.list

ninfo_RNA15_pdb.py VPK_0_34.cg.pdb VPK.hb.list VPK.stack.list a.ninfo
