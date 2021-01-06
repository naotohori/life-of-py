# Convert 1l2x to CG pdb
../../pdb_aa2cg.py ./1l2x.pdb Y ./1l2x.cg.pdb


# Check the CG pdb comparing with the PDB by RNA_cg
bestfit_pdb_chains.py system.pdb 1l2x.cg.pdb

# Output: 
#  0.0006863012572860542
#  [[ 6.91769311e-01 -2.95869811e-01  6.58723216e-01  3.25762153e+02]
#   [ 5.96455058e-01  7.48324065e-01 -2.90262738e-01  3.39790917e+02]
#   [-4.07058454e-01  5.93693648e-01  6.94140668e-01  3.20704309e+02]
#   [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.00000000e+00]]

ninfo_RNA.py --model NHT19 --pdb 1l2x.cg.pdb --hbfile bwyv.hb.dat --tstfile bwyv.tst.dat --end5 P --end3 P --nnhb 3 --exvfile bwyv19.exv.dat bwyv19.ninfo
