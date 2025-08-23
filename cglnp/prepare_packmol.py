#!/usr/bin/env python

import numpy as np
import random as rng

def psf_atom(f, lines, atomid, atomname, resid, resname, chainid):
    if f is None:
        return
    lines.append(f'{atomid:8d} {resname[:4]:4s} {str(resid)[-4:]:4s} {resname[:4]:4s} {atomname[:4]:4s} {str(chainid)[-4:]:4s} {0.0:14.8g} {0.0:13.8g} {0:8d}\n')

#salts = ['Inf']
#replicates = 1
salts = [0.05, 0.1, 0.2, 0.5, 1.0, 'Inf']
replicates = 8

box = 900.0

#molar ratios
molcat = 0.463
moldspc = 0.094
molchol = 0.427
molpeg = 0.016
np_ratio = 6

sigma = 9.98772975127196

log = open(f'./prepare_packmol.log', 'w')

out = open(f'./initial_template.cif', 'w')
out.write(f'data_cg_RNA\n')
out.write('''#
loop_
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
''')

#fpsf = None
fpsf = open(f'./ALC.psf', 'w')
bonds = []     # for PSF
psf_atom_lines = []

rnain = f'./../flucmRNA/CInf/00/last-relax.cg.cif'
catin = './structures/ALC-0315.cif'
dspcin = './structures/dspc.cif'
cholin = './structures/chol.cif'
pegin = './structures/ALC-0159.cif'

with open(rnain, 'r') as f:
    lines = f.readlines()
    for line in lines:
        id = line.split()[0]
        if id.isdigit():
            lsp = line.split()
            x = float(lsp[8])
            y = float(lsp[9])
            z = float(lsp[10])
            out.write(f'{id} {lsp[1]} {lsp[2]} {lsp[3]} {lsp[4]} {lsp[5]} {lsp[6]} {lsp[7]} 0.0 0.0 0.0 {lsp[11]}\n')
            atomname = lsp[1]
            resname = lsp[3]
            rnachainid = int(lsp[5])
            rnaresid = int(lsp[6])
            rnaatomid = int(lsp[0])
            psf_atom(fpsf, psf_atom_lines, rnaatomid, atomname, rnaresid, resname, rnachainid)
            if atomname[0:3] == 'RBB':
                if resname[-1] != '5':
                    bonds.append((rnaatomid-2, rnaatomid))
            else:
                bonds.append((rnaatomid-1, rnaatomid))

chainid = rnachainid
resid = rnaresid
id = rnaatomid

nRNAnt = resid
ncat = int(nRNAnt*np_ratio)
log.write(f'RNA residues: {nRNAnt}\n')
log.write(f'molcat {molcat}\n')
log.write(f'ncat {ncat}\n')
log.flush()
chainid += 1

nbeads_cat = 0
resnames = []
atomnames = []
with open(catin, 'r') as f:
    lines = f.readlines()
    for line in lines:
        atomid = line.split()[0]
        if atomid.isdigit():
            lsp = line.split()
            resname = lsp[1]
            atomname = lsp[-1]
            nbeads_cat += 1
            resnames.append(resname)
            atomnames.append(atomname)

i = 0
while i < ncat:
    resid += 1
    for atom in range(nbeads_cat):
        id += 1
        out.write(f'{id} {resnames[atom]} . ILN {chainid} {chainid} {resid} . {0.0:.3f} {0.0:.3f} {0.0:.3f} {atomnames[atom]}\n')
        psf_atom(fpsf, psf_atom_lines, id, resnames[atom], resid, 'ILN', chainid)
    bonds.append((id-6, id-5))
    bonds.append((id-5, id-4))
    bonds.append((id-4, id-3))
    bonds.append((id-3, id-2))
    bonds.append((id-2, id-1))
    bonds.append((id-1, id))
    i += 1
print(f'Added {i} cationic lipids')
log.write(f'Added {i} cationic lipids\n')

ndspc = int(ncat/molcat*moldspc)
log.write(f'moldspc {moldspc}\n')
log.write(f'ndspc {ndspc}\n')
log.flush()
chainid += 1

nbeads_dspc = 0
resnames = []
atomnames = []
with open(dspcin, 'r') as f:
    lines = f.readlines()
    for line in lines:
        atomid = line.split()[0]
        if atomid.isdigit():
            lsp = line.split()
            resname = lsp[1]
            atomname = lsp[-1]
            nbeads_dspc += 1
            resnames.append(resname)
            atomnames.append(atomname)

i = 0
while i < ndspc:
    resid += 1
    for atom in range(nbeads_dspc):
        id += 1
        out.write(f'{id} {resnames[atom]} . DSPC {chainid} {chainid} {resid} . {0.0:.3f} {0.0:.3f} {0.0:.3f} {atomnames[atom]}\n')
        psf_atom(fpsf, psf_atom_lines, id, resnames[atom], resid, 'DSPC', chainid)
    bonds.append((id-4, id-3))
    bonds.append((id-3, id-4))
    bonds.append((id-3, id-1))
    bonds.append((id-1, id))
    i += 1
print(f'Added {i} dspc lipids')
log.write(f'Added {i} dspc lipids\n')

nchol = int(ncat/molcat*molchol)
log.write(f'molchol {molchol}\n')
log.write(f'nchol {nchol}\n')
log.flush()
chainid += 1

nbeads_chol = 0
resnames = []
atomnames = []
with open(cholin, 'r') as f:
    lines = f.readlines()
    for line in lines:
        atomid = line.split()[0]
        if atomid.isdigit():
            lsp = line.split()
            resname = lsp[1]
            atomname = lsp[-1]
            nbeads_chol += 1
            resnames.append(resname)
            atomnames.append(atomname)

i = 0
while i < nchol:
    resid += 1
    for atom in range(nbeads_chol):
        id += 1
        out.write(f'{id} {resnames[atom]} . CHOL {chainid} {chainid} {resid} . {0.0:.3f} {0.0:.3f} {0.0:.3f} {atomnames[atom]}\n')
        psf_atom(fpsf, psf_atom_lines, id, resnames[atom], resid, 'CHOL', chainid)
    bonds.append((id-2, id-1))
    bonds.append((id-1, id))
    i += 1
print(f'Added {i} chol lipids')
log.write(f'Added {i} chol lipids\n')

npeg = int(ncat/molcat*molpeg)
log.write(f'molpeg {molpeg}\n')
log.write(f'npeg {npeg}\n')
log.flush()

nbeads_peg = 0
resnames = []
atomnames = []
chainnames = []
with open(pegin, 'r') as f:
    lines = f.readlines()
    for line in lines:
        atomid = line.split()[0]
        if atomid.isdigit():
            lsp = line.split()
            resname = lsp[1]
            chainname = lsp[3]
            atomname = lsp[-1]
            nbeads_peg += 1
            resnames.append(resname)
            atomnames.append(atomname)
            chainnames.append(chainname)

i = 0
while i < npeg:
    resid += 1
    chainid += 1
    for atom in range(nbeads_peg):
        id += 1
        if chainnames[atom].startswith('PEG'):
            resid += 1
        out.write(f'{id} {resnames[atom]} . {chainnames[atom]} {chainid} {chainid} {resid} . {0.0:.3f} {0.0:.3f} {0.0:.3f} {atomnames[atom]}\n')
        psf_atom(fpsf, psf_atom_lines, id, resnames[atom], resid, chainnames[atom], chainid)
    bonds.append((id-14, id-13))
    bonds.append((id-14, id-12))
    bonds.append((id-14, id-11))
    bonds.append((id-11, id-10))
    bonds.append((id-10, id-9))
    bonds.append((id-9, id-8))
    bonds.append((id-8, id-7))
    bonds.append((id-7, id-6))
    bonds.append((id-6, id-5))
    bonds.append((id-5, id-4))
    bonds.append((id-4, id-3))
    bonds.append((id-3, id-2))
    bonds.append((id-2, id-1))
    bonds.append((id-1, id))
    i += 1
print(f'Added {i} peg lipids')
log.write(f'Added {i} peg lipids\n')

#Solvate
log.write('Solvating...\n')

n = int((3*(box**3))/(sigma**3))
nsolvent = n - id
vETOH = 0.2
WtETH = (0.997*(1-vETOH)/18.01528)/(0.789*vETOH/46.06844)
Wratio = WtETH*3/10
nethanol = round(nsolvent / (Wratio + 1))
nwater = nsolvent - nethanol
log.write(f'n {n}\n')
log.write(f'nsolvent {nsolvent}\n')
log.write(f'vETOH {vETOH}\n')
log.write(f'WtETH {WtETH}\n')
log.write(f'Wratio {Wratio}\n')
log.write(f'nethanol {nethanol}\n')
log.write(f'nwater {nwater}\n')
log.flush()

chainid += 1
iETOH = 0
for i in range(nethanol):
    resid += 1
    id += 1
    iETOH += 1
    out.write(f'{id} ETH . ETH {chainid} {chainid} {resid} . {0.0:.3f} {0.0:.3f} {0.0:.3f} Th\n')
    psf_atom(fpsf, psf_atom_lines, id, 'ETH', resid, 'ETH', chainid)

chainid += 1
iW = 0
for i in range(nwater):
    resid += 1
    id += 1
    iW += 1
    out.write(f'{id} WAT . WAT {chainid} {chainid} {resid} . {0.0:.3f} {0.0:.3f} {0.0:.3f} W\n')
    psf_atom(fpsf, psf_atom_lines, id, 'WAT', resid, 'WAT', chainid)

log.write(f'Added {iETOH} Ethanols\n')
log.write(f'Added {iW} waters\n')
log.flush()

out.write('''#
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
''')

with open(rnain, 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('covale'):
            out.write(line)

chainid = rnachainid
resid = rnaresid
chainid += 1

for i in range(ncat):
    resid += 1
    out.write(f'covale {chainid} {resid} ILNa {chainid} {resid} ILNb\n')
    out.write(f'covale {chainid} {resid} ILNb {chainid} {resid} ILNc\n')
    out.write(f'covale {chainid} {resid} ILNc {chainid} {resid} ILNd\n')
    out.write(f'covale {chainid} {resid} ILNd {chainid} {resid} ILNe\n')
    out.write(f'covale {chainid} {resid} ILNe {chainid} {resid} ILNf\n')
    out.write(f'covale {chainid} {resid} ILNf {chainid} {resid} ILNg\n')

chainid += 1

for i in range(ndspc):
    resid += 1
    out.write(f'covale {chainid} {resid} DSPC3 {chainid} {resid} DSPC4\n')
    out.write(f'covale {chainid} {resid} DSPC4 {chainid} {resid} DSPC5\n')
    out.write(f'covale {chainid} {resid} DSPC4 {chainid} {resid} DSPC2\n')
    out.write(f'covale {chainid} {resid} DSPC2 {chainid} {resid} DSPC1\n')

chainid += 1

for i in range(nchol):
    resid += 1
    out.write(f'covale {chainid} {resid} CHOL1 {chainid} {resid} CHOL2\n')
    out.write(f'covale {chainid} {resid} CHOL2 {chainid} {resid} CHOL3\n')

for i in range(npeg):
    resid += 1
    chainid += 1
    out.write(f'covale {chainid} {resid} AL1a {chainid} {resid} AL1b\n')
    out.write(f'covale {chainid} {resid} AL1a {chainid} {resid} AL1c\n')
    out.write(f'covale {chainid} {resid} AL1a {chainid} {resid} PEGa\n')
    out.write(f'covale {chainid} {resid} PEGa {chainid} {resid+1} PEG\n')
    for n in range(9):
        resid += 1
        out.write(f'covale {chainid} {resid} PEG {chainid} {resid+1} PEG\n')
    resid += 1
    out.write(f'covale {chainid} {resid} PEG {chainid} {resid+1} PEGT\n')
    resid += 1

out.write('''#''')
out.close()

if fpsf is not None:
    fpsf.write('PSF\n')
    fpsf.write('\n')

    n = len(psf_atom_lines)
    n_solute = 0
    fpsf.write(f'{n:8d} !NATOM\n')
    for l in psf_atom_lines:
        fpsf.write(l)
        if l[9:12] not in ['WAT', 'ETH']:
            n_solute += 1
    fpsf.write('\n')

    n = len(bonds)
    fpsf.write(f'{n:8d} !NBOND: bonds\n')
    for ibond, bond in enumerate(bonds):
        fpsf.write(f'{bond[0]:8d}{bond[1]:8d}')
        if ibond % 4 == 3:
            fpsf.write('\n')
    fpsf.write('\n')

    fpsf.close()


for salt in salts:
    for rep in range(replicates):

        print(f'Working on C{salt} {rep:02d}')

        # Write rna.xyz
        rnain = f'../flucmRNA/C{salt}/{rep:02d}/last-relax.cg.cif'

        nbeads_rna = 0
        with open(rnain, 'r') as f:
            lines = f.readlines()
            for line in lines:
                id = line.split()[0]
                if id.isdigit():
                    nbeads_rna += 1

        rnaxyz = open(f'./C{salt}/{rep:02d}/rna.xyz', 'w')
        rnaxyz.write(f'{nbeads_rna}\n')
        rnaxyz.write(f'flucmRNA\n')
        with open(rnain, 'r') as f:
            lines = f.readlines()
            for line in lines:
                id = line.split()[0]
                if id.isdigit():
                    lsp = line.split()
                    x = float(line.split()[8])
                    y = float(line.split()[9])
                    z = float(line.split()[10])
                    element = line.split()[11]
                    rnaxyz.write(f'{element}  {x} {y} {z}\n')
        rnaxyz.close()

        # Write packmol input file
        pkminp = open(f'./C{salt}/{rep:02d}/packmol.inp', 'w')
        pkminp.write(f'''tolerance {sigma}
filetype xyz
output initial.xyz

# RNA
structure ./rna.xyz
   number 1
   center
   fixed 500. 500. 500. 0. 0. 0.
   radius 0.8
end structure

# Cationic lipid
structure ../../structures/ALC-0315.xyz
   number {ncat}
   pbc 0. 0. 0. 1000. 1000. 1000.
   radius 0.8
end structure

# DSPC
structure ../../structures/dspc.xyz
   number {ndspc}
   pbc 0. 0. 0. 1000. 1000. 1000.
end structure

# Chol
structure ../../structures/chol.xyz
   number {nchol}
   pbc 0. 0. 0. 1000. 1000. 1000.
end structure

# PEG
structure ../../structures/ALC-0159.xyz
   number {npeg}
   pbc 0. 0. 0. 1000. 1000. 1000.
   radius 0.8
end structure

# Ethanol
structure ../../structures/ethanol.xyz
   number {nethanol}
   pbc 0. 0. 0. 1000. 1000. 1000.
   radius 0.5
end structure

# water
structure ../../structures/water.xyz
   number {nwater}
   pbc 0. 0. 0. 1000. 1000. 1000.
   radius 0.5
end structure
''')
        pkminp.close()

        # Write packmol condor file
        condor = open(f'./C{salt}/{rep:02d}/packmol.condor', 'w')
        condor.write(f'''
batch_name = pkmC{salt}_{rep:02d}
FNAME = packmol
ID       = $(Cluster).$(Process)
output   = $(FNAME).$(ID).out
error    = $(FNAME).$(ID).err
log      = $(FNAME).$(Cluster).log

universe = vanilla

getenv = True
executable    = /users/paznh/local/bin/packmol
input         = ./packmol.inp

request_memory = 4 GB
Request_CPUs = 1
Request_GPUs = 0

queue
''')
        condor.close()


