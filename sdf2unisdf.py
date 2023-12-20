#!/usr/bin/env python
# coding: utf-8

from rdkit import Chem
from rdkit.Chem import build_torsion_tree
import os,sys,string,argparse
from optparse import OptionParser

parser = argparse.ArgumentParser(description="Prepare Uni-Dock SDF.\n")
parser.add_argument('input',metavar='<input>',help="input sdf file")
parser.add_argument('output',metavar='<output>',help="output sdf file")
args = parser.parse_args()
ifile = args.input
ofile = args.output

suppl = Chem.SDMolSupplier(ifile, removeHs=False)
ofile = Chem.SDWriter(ofile)

for mol in suppl:
   if mol is not None:
       runner = build_torsion_tree.TopologyBuilder(mol)
       runner.build_molecular_graph()
       fragment_info_string, torsion_info_string, atom_info_string = runner.get_sdf_torsion_tree_info()
       mol.SetProp("fragInfo", fragment_info_string)
       mol.SetProp("torsionInfo", torsion_info_string)
       mol.SetProp("atomInfo", atom_info_string)
       ofile.write(mol)
ofile.close()
