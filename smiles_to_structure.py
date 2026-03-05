#!/usr/bin/env python

from openbabel import openbabel as ob

def smiles_to_structure(smiles, pH, name):
    conv = ob.OBConversion()
    conv.SetInFormat("smi")

    mol = ob.OBMol()
    conv.ReadString(mol, smiles)

    ## Below, AddHydrogens then CorrectForPH does not work.
    ## Add hydrogens
    #mol.AddHydrogens()
    # Adjust protonation states depending on pH
    #mol.CorrectForPH(pH)


    # Instead, do all by AddHydrogens, providing three arguments.
    # 1st: polaronly
    # 2nd: correctForPH
    # 3rd: pH
    mol.AddHydrogens(False, True, pH)


    # Build 3D structure
    builder = ob.OBBuilder()
    builder.Build(mol)

    # Optimisation
    #ff = ob.OBForceField.FindForceField("MMFF94")
    #ff = ob.OBForceField.FindForceField("GAFF")
    ff = ob.OBForceField.FindForceField("AM")
    ff.Setup(mol)
    ff.ConjugateGradients(500)
    ff.GetCoordinates(mol)

    # smiles output
    conv.SetOutFormat("smi")
    conv.WriteFile(mol, name + ".smiles.txt")

    # mol2 output
    conv.SetOutFormat("mol2")
    conv.WriteFile(mol, name + ".mol2")

    # PDB output
    conv.SetOutFormat("pdb")
    conv.WriteFile(mol, name + ".pdb")

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Process input arguments")

    parser.add_argument("smiles", type=str, help="SMILES string")
    parser.add_argument("--name", type=str, default="new", help='Output file prefix (default: "new")')
    parser.add_argument("--pH", type=float, default=7.4, help="pH value (default: 7.4)")

    args = parser.parse_args()

    smiles_to_structure(args.smiles, args.pH, args.name)

