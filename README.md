# protein-contact-map-generator
A python based, easy to use contact map generator.  Comments and suggestions are welcome.

# Usage：
  python PCMG.py -pdb_id 1B0B -display
  
  or python PCMG.py -pdb_path c:\1B0B.pdb -seq_path c:\1B0B.fastav -o c:\1B0B.cm
## optional arguments:
  
  -pdb_id     pdb id to download pdb file and sequence file 
  
  -pdb_path   pdb file used to caculate distance between residues
  
  -seq_path   sequence file used to renumber residues consecutively
  
  -chain_id   chain to generate
  
  -display    show the contact map by pickle
  
  -o          store the contact map using pickle
  
  

  

  
