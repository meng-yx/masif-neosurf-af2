#!/usr/bin/python
import Bio
from Bio.PDB import * 
import sys
import importlib
import os
import urllib.request

from default_config.masif_opts import masif_opts
# Local includes
from input_output.protonate import protonate

if len(sys.argv) <= 1: 
    print("Usage: "+sys.argv[0]+" PDBID_A_B")
    print("A or B are the chains to include in this pdb.")
    sys.exit(1)

if not os.path.exists(masif_opts['raw_pdb_dir']):
    os.makedirs(masif_opts['raw_pdb_dir'])

if not os.path.exists(masif_opts['tmp_dir']):
    os.mkdir(masif_opts['tmp_dir'])

in_fields = sys.argv[1].split('_')
pdb_id = in_fields[0]

# Download pdb - use direct HTTPS download instead of BioPython PDBList
# which has issues with the FTP server
import ssl

pdb_filename = os.path.join(masif_opts['tmp_dir'], 'pdb{}.ent'.format(pdb_id.lower()))
pdb_url = 'http://files.rcsb.org/download/{}.pdb'.format(pdb_id.upper())
print('Downloading PDB structure from {}...'.format(pdb_url))

# Create SSL context
ctx = ssl.create_default_context()

try:
    with urllib.request.urlopen(pdb_url, context=ctx) as response:
        with open(pdb_filename, "wb") as f:
            f.write(response.read())
    print('PDB file downloaded successfully')
except Exception as e:
    print('Error downloading PDB file: {}'.format(e))
    sys.exit(1)


##### Protonate with reduce, if hydrogens included.
# - Always protonate as this is useful for charges. If necessary ignore hydrogens later.
protonated_file = masif_opts['raw_pdb_dir']+"/"+pdb_id+".pdb"
protonate(pdb_filename, protonated_file)
pdb_filename = protonated_file

