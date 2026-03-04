#!/usr/bin/env python

from openbabel import openbabel as ob
import sys

if len(sys.argv) not in (2,3):
    print("Usage: SCRIPT <SMILES string> [pH value]")
    print("       If pH is not given, default pH = 7.4 is used.")
    sys.exit(0)

s = sys.argv[1]

if len(sys.argv) == 3:
    pH = float(sys.argv[2])
else:
    pH = 7.4

obc = ob.OBConversion()

obc.SetInAndOutFormats('smi', 'smi')

def get_smi_with_pH(smi, pH=pH):
    obmol = ob.OBMol()
    obc.ReadString(obmol, smi)
    obmol.CorrectForPH(pH)
    return obc.WriteString(obmol)

print(get_smi_with_pH(s))
