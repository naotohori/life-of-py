#!/usr/bin/env python

import sys
import re


if len(sys.argv) != 3:
    print("Usage: SCRIPT [DSSR file] [Output prefix]")
    print("This script read the DSSR file, and generate following files:")
    print("    PREFIX.dssr.out .... Summarize hbond information from DSSR")
    print("    PREFIX.hb.native.dat .... hbond file to semiexplicit model, exclusively native hbonds from DSSR" )
    print("    PREFIX.hbond.dat .... hbond file to semiexplicit model, including additional all possible non-native base pairs" )
    sys.exit(2)

file_dssr = sys.argv[1]
file_out_prefix = sys.argv[2]


""""""""""""""""""""""""""""""""""""""
"""""""""     Read DSSR      """""""""
""""""""""""""""""""""""""""""""""""""

flg_read = False
count_arrow = 0
seq = ''

for l in open(file_dssr):

    if not flg_read:
        lsp = l.split()

        if len(lsp) < 5:
            continue

        if lsp[0] != 'Secondary' or lsp[1] != 'structures' or lsp[2] != 'in' or lsp[3] !='dot-bracket' or lsp[4] != 'notation':
            pass
        else:
            flg_read = True

        continue

    if l[0] in ('.', '(', ')', '[', ']'):
        break

    if l[0] == '>':
        count_arrow += 1
        if count_arrow == 1:
            continue
        else:
            break

    seq += l.strip()

print('Sequence: %s' % seq)
print('#nt: %i' % len(seq))


flg_read = False
lines_base_pair = []

for l in open(file_dssr):

    if not flg_read:
        lsp = l.split()

        if len(lsp) < 5:
            continue

        if lsp[0] != 'List' or lsp[1] != 'of' or lsp[3] != 'base' or lsp[4] !='pairs':
            pass
        else:
            nbp = int(lsp[2])
            flg_read = True

        continue

    if len(l.strip()) == 0:
        break

    lines_base_pair.append(l)

# Remove the first line 
# "    nt1            nt2            bp  name        Saenger   LW   DSSR"
lines_base_pair.pop(0)

# Check the number of lines
if len(lines_base_pair) % 6 != 0:
    print("Error: len(lines_base_pair) % 6 != 0")
    sys.exit(2)

if nbp != len(lines_base_pair) // 6:
    print("Error: nbp != len(lines_base_par) // 6")
    sys.exit(2)

print('Number of base pairs: %i' % nbp)


""""""""""""""""""""""""""""""""""""""
"""""""""     Output         """""""""
""""""""""""""""""""""""""""""""""""""

f_out = open(file_out_prefix+'.dssr.out','w')
f_hbnative = open(file_out_prefix+'.hb.native.dat','w')
f_hb = open(file_out_prefix+'.hbond.dat','w')

def write_both(s):
    f_hbnative.write(s)
    f_hb.write(s)

native_pairs = []

for ibp in range(nbp):

    """ Extract information """

    l_main = lines_base_pair[6*ibp]
    '''   1 U3             G34            U-G --          n/a       cHW  cM-W'''
    l_hbonds = lines_base_pair[6*ibp+3]
    '''H-bonds[2]: "O4(carbonyl)-N1(imino)[2.80],O4(carbonyl)-N2(amino)[2.91]"'''
    #print(l_main)
    #print(l_hbonds)

    l_main_sp = l_main.split()
    if ibp != int(l_main_sp[0])-1:
        print('Error: ibp != int(l_main_sp[0])')
        sys.exit(2)

    nt1 = l_main_sp[1]
    nt2 = l_main_sp[2]
    bpinfo_bp = l_main_sp[3]
    bpinfo_name = l_main_sp[4]
    bpinfo_Saenger = l_main_sp[5]
    bpinfo_LW = l_main_sp[6]
    bpinfo_DSSR = l_main_sp[7]

    if len(bpinfo_bp) > 3:
        print("Warning: len(bpinfo_bp) > 3")
    if len(bpinfo_name) > 8:
        if bpinfo_name != "~rHoogsteen":
            print("Warning: len(bpinfo_name) > 8")
        ''' 8 "~Sheared", "Platform" '''
        ''' 11 "~rHoogsteen" '''
    if len(bpinfo_Saenger) > 9:
        print("Warning: len(bpinfo_Saenger) > 9")
    if len(bpinfo_LW) > 3:
        print("Warning: len(bpinfo_LW) > 3")
    if len(bpinfo_DSSR) > 4:
        print("Warning: len(bpinfo_DSSR) > 4")


    """ Write specific to out file """

    f_out.write('%5i %5s %5s' % (ibp+1, nt1, nt2))

    r = re.search(r'H-bonds\[(\d+)\]: "(.*)"', l_hbonds)
    nHB = int(r.group(1))
    details = r.group(2)
    #print(nHB)
    #print(details)

    f_out.write('  %i  ' % nHB)
    f_out.write(' %3s' % bpinfo_bp)
    f_out.write(' %8s' % bpinfo_name)
    f_out.write(' %9s' % bpinfo_Saenger)
    f_out.write(' %3s' % bpinfo_LW)
    f_out.write(' %4s' % bpinfo_DSSR)



    """ Write specific to native.hb.dat file """

    write_both('NAT')
    if bpinfo_name in ('WC', 'Wobble'):
        write_both('%s-%s' % (bpinfo_bp[0], bpinfo_bp[2]))
    else:
        write_both('IVE')

    r = re.fullmatch(r'\D+(\d+)', nt1)
    if r is None:
        print('Error: cannot recognize nt1, %s' % nt1)
    nnt1 = int(r.group(1))
    write_both(' %5i' % nnt1)

    r = re.fullmatch(r'\D+(\d+)', nt2)
    if r is None:
        print('Error: cannot recognize nt1, %s' % nt2)
    nnt2 = int(r.group(1))
    write_both(' %5i' % nnt2)

    write_both('  ')

    native_pairs.append( (nnt1, nnt2) )
    
    """ Write common field (ATOMS) to the both files """

    for detail in details.split(','):
        '''
        Note:
        The connector can be either "-" or "*"
        e.g.,
        H-bonds[4]: "N2(amino)-O4(carbonyl)[3.33],N3*O4(carbonyl)[3.23],O2'(hydroxyl)-O5'[2.84],O2'(hydroxyl)-O4'[3.12]"
        '''
        r = re.search(r'([^\(\[]*)[^-]*[-\*]([^\(\[]*)', detail)
        f_out.write(' %s %s' % (r.group(1), r.group(2)))
        write_both(' %s %s' % (r.group(1), r.group(2)))

    f_out.write(' TER\n')
    write_both(' TER\n')


f_out.close()
f_hbnative.close()

###############################################################################

NNHB_NT_SEP = 5

atoms_AU = 'N6 O4 N1 N3 TER'
atoms_GC = 'N1 N3 N2 O2 O6 N4 TER'
atoms_GU = 'N1 O2 O6 N3 TER'

nseq = len(seq)
print('length = %i' % nseq)


for i in range(2, nseq):
    nt_i = i + 1

    for j in range(i+NNHB_NT_SEP, nseq):
        nt_j = j + 1

        if (nt_i, nt_j) in native_pairs:
            continue

        if seq[i] == 'G' and seq[j] == 'C':
            f_hb.write('G-C    %4i %4i  %s\n' % (nt_i, nt_j, atoms_GC))
        if seq[i] == 'C' and seq[j] == 'G':
            f_hb.write('G-C    %4i %4i  %s\n' % (nt_j, nt_i, atoms_GC))
        if seq[i] == 'A' and seq[j] == 'U':
            f_hb.write('A-U    %4i %4i  %s\n' % (nt_i, nt_j, atoms_AU))
        if seq[i] == 'U' and seq[j] == 'A':
            f_hb.write('A-U    %4i %4i  %s\n' % (nt_j, nt_i, atoms_AU))
        if seq[i] == 'G' and seq[j] == 'U':
            f_hb.write('G-U    %4i %4i  %s\n' % (nt_i, nt_j, atoms_GU))
        if seq[i] == 'U' and seq[j] == 'G':
            f_hb.write('G-U    %4i %4i  %s\n' % (nt_j, nt_i, atoms_GU))

f_hb.close()
