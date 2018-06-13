import Bio.PDB
import numpy as np
import requests
from itertools import islice
import os
import warnings
import sys
import argparse 

DIC = {'ALA': 'A', 'LEU': 'L', 'ARG': 'R', 'LYS': 'K', 'ASN': 'N', 'MET': 'M', 'ASP': 'D', 'PHE': 'F', 'CYS': 'C',
       'PRO': 'P', 'GLN': 'Q', 'SER': 'S', 'GLU': 'E', 'THR': 'T', 'GLY': 'G', 'TRP': 'W', 'HIS': 'H', 'TYR': 'Y',
       'ILE': 'I', 'VAL': 'V','UNK':'X'}
def calc_gap(chain,seq):
    pdb_seq = []
    seqs = chain.child_list[0].get_id()[1]
    gap = 0
    for i in chain:
        while len(pdb_seq) != i.get_id()[1] - seqs:
            pdb_seq.append('X')
        pdb_seq.append(DIC[i.get_resname()])

    for i, j in enumerate(seq):
        match = True
        for k, l in enumerate(pdb_seq):
            if l != 'X' and seq[k + i] != 'X' and seq[k + i] != l:
                match = False
                break
        if match:
            gap = i - seqs + 1
            break
    else:
        raise Exception("inconsistent between pdb files and sequence files")
    return gap

def clear_pdb(infile,outfile):
    PDBtxt = ''
    with open(infile, 'rU') as f: 
        for line in f.read().splitlines():
            if line.startswith("END"):
                line = line.replace("ENDMDL", "END   ")
                PDBtxt += line + '\n'
                break
            if line.startswith("ATOM  ") or line.startswith("TER"):
                PDBtxt += line + '\n'
    with open(outfile, 'w') as f:
        f.write(PDBtxt)


def calc_residue_dist(residue_one, residue_two):
    c1 = 'CB' if residue_one.get_resname() != 'GLY' else 'CA'
    c2 = 'CB' if residue_two.get_resname() != 'GLY' else 'CA'
  
    if c1 not in residue_one:
        c1=residue_one.child_list[0].id
    if c2 not in residue_two:
        c2=residue_two.child_list[0].id

    diff_vector = residue_one[c1].coord - residue_two[c2].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

def generate_contact_map(pdb_path,seq_path):
    if not os.path.exists(pdb_path) :
        raise Exception("No such pdb file :{0}".format(pdb_path))
    elif not os.path.exists(seq_path):
        raise Exception("No such sequence file :{0}".format(seq_path))
    _,pdb_fname=os.path.split(pdb_path)   
    _,pdb_name=os.path.splitext(pdb_path)   
    with open(seq_path) as f:
        contents = f.read()
    sequence=''.join(contents.split('\n')[1:-1])
    length=len(sequence)
    
    clear_pdb(pdb_path,'cache/'+pdb_fname)
    chain = Bio.PDB.PDBParser().get_structure(pdb_name,'cache/'+pdb_fname)[0]['A'] #
    gap=calc_gap(chain,sequence)
    contact_map = np.zeros((length, length))
    for row in range(length):
        for col in range(length):
            if chain.has_id((' ', row + 1-gap, ' ')) and chain.has_id((' ', col + 1-gap, ' ')):
                contact_map[row, col] = calc_residue_dist(chain.__getitem__(row + 1-gap), chain.__getitem__(col + 1-gap))
            else:
                contact_map[row, col] = -1
    return contact_map
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument("-pdb_id")
    parser.add_argument("-pdb_path")
    parser.add_argument("-seq_path")
    parser.add_argument("-display", action="store_true")
    parser.add_argument("-o")
    
    args = parser.parse_args()

    if args.pdb_path and args.seq_path:
        cm=generate_contact_map(args.pdb_path,args.seq_path) #E:\s2\sequence\1B0BA.seq
        
    if args.display:
        import matplotlib.pyplot as plt
        plt.imshow(cm,cmap='gray_r')
        plt.show()
    if args.o:
        if os.path.split(args.o)[0]!='' and  not os.path.exists(os.path.split(args.o)[0]):
            raise Exception('path don\'t exists')
        import pickle
        with open(args.o,'wb') as f:
            pickle.dump(cm,f)
        


    pass