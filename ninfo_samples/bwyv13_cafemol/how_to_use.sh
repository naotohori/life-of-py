## ninfo for a single stranded sequence without base pairs.
#  e.g. polyA
ninfo_RNA_cafemol.py --seq "AAAAAAAAAA" --model DT13 A10.ninfo


## ninfo with specific base pairs (hbfile). e.g., BWYV pseudoknot 

# 1. Prepare a list of hydrogen bonds (e.g. bwyv.hb.list)
#    Thiss can be generated usingh DSSP (with manual editing).
#   e.g.) x3dna-dssr -i=BWYV_0_29.pdb -o=BWYV.dssr.out --prefix=BWYV.dssr --non-pair

# 2. Generate cgpdb.
pdb_aa2cg.py BWYV_0_29.pdb n BWYV_0_29.cg.pdb

# 3. Run ninfo_RNA_cafemol.py to generate an ninfo file.
ninfo_RNA_cafemol.py --model DT13 --pdb BWYV_0_29.cg.pdb --hbfile bwyv.hb.list bwyv_RNA13.ninfo

