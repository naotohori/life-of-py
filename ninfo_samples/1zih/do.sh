mkdir dssr
x3dna-dssr -i=1zih_1.pdb -o=dssr/1zih.dssr.out --prefix=dssr/1zih.dssr --non-pair

#pdb_cg2pdb.py 1zih_1.pdb 1zih.cg.pdb
pdb_aa2cg.py 1zih_1.pdb 1zih.cg.pdb

# Edit 1zih.cg.pdb to add residue 0 and 13

ninfo_RNA13_pdb.py --pdb 1zih.cg.pdb --hbfile 1zih.hb.list 1zih.ninfo
