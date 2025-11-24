#!/usr/bin/env python3
"""
CLI script to download AF2 equivalent regions of PDB structures.

For a given PDB structure, this script:
- For each chain:
    - Queries PDB API for uniprot_id for this chain
    - Queries AF2 database for the AF2 model corresponding to this chain, then downloads it
    - Loads AF2 model and PDB structures into python
    - Gets the sequence of PDB structure as well as AF2 model, computes sequence alignment
    - Deletes residues in AF2 model that do not align with PDB structure
    - Aligns the AF2 model to the PDB structure to reconstruct the complex
"""

import argparse
import os
import sys
import json
import requests
import tempfile
import subprocess
import re
import copy
from pathlib import Path
from textwrap import wrap
from io import StringIO

from Bio.PDB import PDBList, PDBParser, PDBIO, Structure, Model
from Bio.PDB import Superimposer
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio import pairwise2
from Bio.SeqUtils import seq1


# ============================================================================
# Helper Functions
# ============================================================================

def get_uniprot_id_from_pdb_chain(pdb_id, chain_id, verbose=False):
    """
    Queries the PDBe API to get the UniProt accession for a given PDB ID and chain.

    Args:
        pdb_id (str): The PDB ID (e.g., '7MON').
        chain_id (str): The chain identifier (e.g., 'A').
        verbose (bool): Whether to print verbose messages.

    Returns:
        str or None: The UniProt ID if found, else None.
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    response = requests.get(url)
    if response.status_code != 200:
        if verbose:
            print(f"Failed to fetch data for PDB {pdb_id}")
        return None
    data = response.json()
    mappings = data.get(pdb_id.lower(), {}).get("UniProt", {})
    for uniprot_id, mapping in mappings.items():
        for segment in mapping.get("mappings", []):
            if segment.get("chain_id") == chain_id:
                return uniprot_id
    if verbose:
        print(f"UniProt ID not found for chain {chain_id} in PDB {pdb_id}")
    return None


def fetch_af2_model(uniprot_id, verbose=False):
    """
    Downloads an AlphaFold2 (AF2) model by querying the AlphaFold DB API using the UniProt ID,
    saves the corresponding PDB file to a temporary location, and loads it into Python using Biopython.

    Args:
        uniprot_id (str): The UniProt ID (e.g., 'Q8NB16').
        verbose (bool): Whether to print verbose messages.

    Returns:
        tuple: (pdb_url, structure)
            pdb_url (str or None): The URL to the downloaded PDB file or None if not found.
            structure: The Biopython Structure object parsed from the downloaded AF2 PDB file,
                       or None if download/parsing failed.
    """
    # Query the AlphaFold API to get pdbUrl
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        api_response = requests.get(api_url)
        if api_response.status_code != 200:
            if verbose:
                print(f"Failed to query AlphaFold API for {uniprot_id} (HTTP {api_response.status_code})")
            return None, None
        api_json = api_response.json()
        if not api_json or 'pdbUrl' not in api_json[0]:
            if verbose:
                print(f"No pdbUrl found in AlphaFold API response for {uniprot_id}")
            return None, None
        pdb_url = api_json[0]['pdbUrl']
    except Exception as e:
        if verbose:
            print(f"Error fetching/parsing AlphaFold API response: {e}")
        return None, None

    # Download the PDB file to a temporary location
    try:
        pdb_response = requests.get(pdb_url)
        if pdb_response.status_code != 200:
            if verbose:
                print(f"Failed to download PDB file for {uniprot_id} at {pdb_url} (HTTP {pdb_response.status_code})")
            return pdb_url, None
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=True) as tmp_fp:
            tmp_fp.write(pdb_response.content)
            tmp_fp.flush()
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(uniprot_id, tmp_fp.name)
            return pdb_url, structure
    except Exception as e:
        if verbose:
            print(f"Failed to download or parse PDB file for {uniprot_id}: {e}")
        return pdb_url, None


def get_protein_sequence(structure):
    """
    Extracts the protein sequence from a Biopython structure.

    Args:
        structure: Biopython Structure object

    Returns:
        dict: A mapping from chain ID to its one-letter protein sequence.
    """
    # Build map from 3-letter (upper) codes to 1-letter codes
    _3to1 = {k.upper(): v for k, v in protein_letters_3to1.items()}

    sequences = {}
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                if is_aa(residue, standard=True):
                    resname = residue.get_resname().upper()
                    try:
                        seq += _3to1[resname]
                    except KeyError:
                        # Skip residues that can't be translated (e.g., unknowns)
                        continue
            if seq:
                sequences[chain.id] = seq
        break  # Only first model (usually sufficient for X-ray/NMR)
    return sequences


def get_single_chain_structure(structure, chain_id):
    """
    Extract a single chain from a structure.

    Args:
        structure: Biopython Structure object
        chain_id: Chain ID to extract

    Returns:
        Biopython Structure object containing only the selected chain
    """
    builder = StructureBuilder()
    builder.init_structure("filtered")
    builder.init_model(0)
    # Add only the selected chain to the new structure
    builder.structure[0].add(structure[0][chain_id])
    return builder.get_structure()


def run_needle_alignment(
    test_seq: str,
    ref_seq: str,
    test_id: str = "test_seq",
    ref_id: str = "ref_seq",
    gapopen: float = 10.0,
    gapextend: float = 0.5,
    needle_exe: str = "needle",
    aformat: str = "pair",
) -> str:
    """
    Run EMBOSS needle to align two protein sequences.

    Parameters
    ----------
    test_seq, ref_seq : str
        Protein sequences (plain one-letter codes, no FASTA header).
    test_id, ref_id : str
        Sequence identifiers for FASTA headers.
    gapopen : float
        Gap opening penalty for needle.
    gapextend : float
        Gap extension penalty for needle.
    needle_exe : str
        Name or path to the 'needle' executable.
    aformat : str
        Output alignment format (e.g. 'pair', 'markx3', 'markx10', 'fasta').

    Returns
    -------
    str
        The alignment output from needle (text).
    """

    def to_fasta(seq: str, header: str) -> str:
        # Wrap sequence at 60 columns for nice FASTA
        wrapped = "\n".join(wrap(seq.replace("\n", "").strip(), 60))
        return f">{header}\n{wrapped}\n"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        a_file = tmpdir / "a_seq.fasta"
        b_file = tmpdir / "b_seq.fasta"

        # Write FASTA files
        a_file.write_text(to_fasta(test_seq, test_id))
        b_file.write_text(to_fasta(ref_seq, ref_id))

        # Ask needle to write to stdout
        cmd = [
            needle_exe,
            "-asequence", str(a_file),
            "-bsequence", str(b_file),
            "-gapopen", str(gapopen),
            "-gapextend", str(gapextend),
            "-aformat", aformat,
            "-stdout",           # send alignment to stdout
            "-auto",             # no interactive prompts
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,  # we'll handle errors ourselves
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"needle failed with code {result.returncode}.\n"
                f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
            )

        return result.stdout


def get_match_seq_range(
    alignment_text: str,
    test_id: str = "test_seq",
    ref_id: str = "ref_seq",
):
    """
    From a needle 'pair' alignment, return the residue indices in ref_seq
    that align to the first and last *aligned* residues in test_seq.

    "Aligned residue" here means a column where BOTH test and ref have
    a non-gap amino acid (not '-').

    Returns
    -------
    (int | None, int | None)
        (first_ref_idx, last_ref_idx), **0-based** indices in ref_seq.
        If no aligned residue pair exists, returns (None, None).
    """

    aligned_test = []
    aligned_ref = []

    # Collect the aligned strings for test and ref
    for line in alignment_text.splitlines():
        line = line.rstrip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        label = parts[0] if parts else ""

        if label == test_id or label == ref_id:
            # Grab the contiguous alignment chunk: letters + '-'
            m = re.search(r'([A-Z\-]+)', line)
            if not m:
                continue
            seq_chunk = m.group(1)

            if label == test_id:
                aligned_test.append(seq_chunk)
            elif label == ref_id:
                aligned_ref.append(seq_chunk)

    aligned_test = "".join(aligned_test)
    aligned_ref = "".join(aligned_ref)

    if len(aligned_test) != len(aligned_ref):
        raise ValueError("Aligned test and ref sequences have different lengths.")

    # i_ref will be 0-based index in the ungapped ref sequence
    i_test = -1  # 0-based index in ungapped test seq, start at -1 before first residue
    i_ref = -1   # 0-based index in ungapped ref seq

    first_ref_idx = None
    last_ref_idx = None

    for a_t, a_r in zip(aligned_test, aligned_ref):
        if a_t != "-":
            i_test += 1
        if a_r != "-":
            i_ref += 1

        # A true aligned residue pair: both non-gaps
        if a_t != "-" and a_r != "-":
            if first_ref_idx is None:
                first_ref_idx = i_ref  # already 0-based
            last_ref_idx = i_ref       # update as we go

    return first_ref_idx, last_ref_idx


def filter_residue_range(structure, start, end):
    """
    Filter the given Bio.PDB structure so that only residues with resseq in [start, end] and hetflag == ' ' remain.
    This function modifies the structure in-place and returns it for convenience.
    """

    for model in structure:
        for chain in model:
            all_res = list(chain)
            for residue in all_res:
                hetflag, resseq, icode = residue.id
                # Only keep ' ' (blank) hetflag (standard residues) and in range
                if hetflag != " ":
                    chain.detach_child(residue.id)
                elif not (start <= resseq <= end):
                    chain.detach_child(residue.id)

    return structure


def align_af2_to_pdb(chain_pdb_struct, chain_af2_struct_filtered):
    """
    Align the filtered AF2 model to the PDB structure using Biopython's Superimposer,
    by performing sequence alignment and only using aligned residues for superimposition.

    Parameters:
        chain_pdb_struct: Bio.PDB.Structure.Structure - structure containing the PDB chain
        chain_af2_struct_filtered: Bio.PDB.Structure.Structure - structure containing the filtered AF2 chain

    Returns:
        chain_af2_struct_aligned: Bio.PDB.Structure.Structure - AF2 structure aligned to the PDB structure
        rmsd: float - RMSD value after alignment (in Angstroms)
    """

    # Deepcopy the AF2 structure so the original is not rotated in-place
    af2_struct_aligned = copy.deepcopy(chain_af2_struct_filtered)

    # Get chains for alignment (use first chain in each structure)
    pdb_chain = next(chain for chain in chain_pdb_struct[0])
    af2_chain = next(chain for chain in af2_struct_aligned[0])

    # Convert residue objects to lists and sequences for the relevant chains
    pdb_residues = [res for res in pdb_chain if res.id[0] == ' ' and 'CA' in res]
    af2_residues = [res for res in af2_chain if res.id[0] == ' ' and 'CA' in res]

    # Derive sequences from these residues
    def get_seq_from_residues(residues):
        return "".join(seq1(res.get_resname()) for res in residues)

    pdb_seq = get_seq_from_residues(pdb_residues)
    af2_seq = get_seq_from_residues(af2_residues)

    # Run global sequence alignment
    aln = pairwise2.align.globalxx(pdb_seq, af2_seq, one_alignment_only=True)[0]
    aln_pdb = aln.seqA
    aln_af2 = aln.seqB

    # Select only aligned (non-gap) C-alpha atoms
    pdb_ca_atoms = []
    af2_ca_atoms = []
    pdb_idx = 0
    af2_idx = 0

    for a, b in zip(aln_pdb, aln_af2):
        if a != '-' and b != '-':
            pdb_ca_atoms.append(pdb_residues[pdb_idx]['CA'])
            af2_ca_atoms.append(af2_residues[af2_idx]['CA'])
        if a != '-':
            pdb_idx += 1
        if b != '-':
            af2_idx += 1

    # Ensure we have matching pairs
    if len(pdb_ca_atoms) != len(af2_ca_atoms):
        raise ValueError(f"After alignment, found {len(pdb_ca_atoms)} matched C-alpha pairs but expected lengths to agree.")

    # Superimpose AF2 (moving) onto PDB (fixed)
    sup = Superimposer()
    sup.set_atoms(pdb_ca_atoms, af2_ca_atoms)
    sup.apply(af2_chain.get_atoms())

    rmsd = sup.rms

    return af2_struct_aligned, rmsd


def process_af2_alignment_single_chain(structure, pdb_id, chain, verbose=False):
    """
    Download the equivalent AF2 model of a chain, align it to the PDB structure,
    and return relevant alignment objects/results.

    Args:
        structure: Biopython structure object for the PDB.
        pdb_id: PDB ID string (e.g., '7MON').
        chain: Chain ID (str, e.g., 'A').
        verbose (bool): Whether to print verbose messages.

    Returns:
        (chain_id, aligned_af2_chain_struct, rmsd, af2_pdb_url)
    """
    # Get the PDB and AF2 structure for a given chain
    chain_pdb_struct = get_single_chain_structure(structure, chain)
    chain_uniprot_id = get_uniprot_id_from_pdb_chain(pdb_id, chain, verbose=verbose)
    chain_af2_pdb_url, chain_af2_struct = fetch_af2_model(chain_uniprot_id, verbose=verbose)

    # Get the PDB and AF2 sequence for a given chain
    chain_pdb_seq = get_protein_sequence(chain_pdb_struct)
    chain_af2_seq = get_protein_sequence(chain_af2_struct)

    # Get the residue range in the AF2 model that matches the PDB sequence
    alignment_text = run_needle_alignment(chain_pdb_seq[chain], chain_af2_seq['A'], test_id="pdb", ref_id="af2")
    match_range = get_match_seq_range(alignment_text, test_id="pdb", ref_id="af2")

    # Filter the AF2 model to only include the matched residues
    chain_af2_struct_filtered = filter_residue_range(chain_af2_struct, match_range[0], match_range[1])

    # Align the filtered AF2 model to the PDB structure
    chain_af2_struct_aligned, af2_to_pdb_rmsd = align_af2_to_pdb(chain_pdb_struct, chain_af2_struct_filtered)

    return chain, chain_af2_struct_aligned, af2_to_pdb_rmsd, chain_af2_pdb_url


def process_af2_alignment(structure, pdb_id, chain_ids, verbose=False):
    """
    For all chain IDs in chain_ids, perform AF2 alignment, and collect results.

    Args:
        structure: Biopython structure object for the PDB.
        pdb_id: PDB ID string (e.g., '7MON').
        chain_ids: List of chain IDs (e.g., ['A','B']).
        verbose (bool): Whether to print verbose messages.

    Returns:
        full_structure: Structure object containing all aligned AF2 chains together as separate chains.
        results: Dictionary mapping chain ID to a dictionary with RMSD value and AF2 PDB URL.
    """
    # Prepare a new structure object to hold all aligned chains
    full_structure = Structure.Structure('aligned_af2')
    model = Model.Model(0)
    full_structure.add(model)

    results = {}

    for chain in chain_ids:
        chain_id, aligned_chain_struct, rmsd, chain_af2_pdb_url = process_af2_alignment_single_chain(structure, pdb_id, chain, verbose=verbose)
        # We expect aligned_chain_struct to have a single chain (id, etc.)
        # Add this chain to our model
        for ch in aligned_chain_struct.get_chains():
            # Remove existing chains in the model with same id to avoid duplicates
            if chain_id in model:
                model.detach_child(chain_id)
            ch.id = chain_id  # Make sure the chain id matches requested chain id
            model.add(ch)
        results[chain_id] = {
            "rmsd": rmsd,
            "af2_pdb_url": chain_af2_pdb_url
        }

    return full_structure, results


# ============================================================================
# Main Workflow
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Download and align AF2 models to PDB structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python PDB_to_AF2.py 7MON_A_B --pdb-out-dir ./data/pdb --rmsd-out-dir ./data/rmsd

This will:
  - Download PDB 7MON
  - Process chain A as interactor 1 and chain B as interactor 2
  - Save aligned AF2 models to ./data/pdb/7MON_A.pdb and ./data/pdb/7MON_B.pdb
  - Save RMSD results to ./data/rmsd/7MON_rmsd.json
        """
    )
    
    parser.add_argument(
        'pdb_complex_id',
        type=str,
        help='PDB complex ID in format PDBID_CHAIN1_CHAIN2 (e.g., 7MON_A_B)'
    )
    
    parser.add_argument(
        '--pdb-out-dir',
        type=str,
        default='./data/pdb',
        help='Output directory for PDB files (default: ./data/pdb)'
    )
    
    parser.add_argument(
        '--rmsd-out-dir',
        type=str,
        default='./data/rmsd',
        help='Output directory for RMSD JSON files (default: ./data/rmsd)'
    )
    
    parser.add_argument(
        '--temp-dir',
        type=str,
        default='.',
        help='Temporary directory for downloading PDB files (default: current directory)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output (default: False)'
    )
    
    args = parser.parse_args()
    
    # Parse PDB complex ID
    parts = args.pdb_complex_id.split('_')
    if len(parts) != 3:
        print(f"Error: PDB complex ID must be in format PDBID_CHAIN1_CHAIN2 (e.g., 7MON_A_B)")
        sys.exit(1)
    
    pdb_id = parts[0]
    interactor_1 = parts[1]
    interactor_2 = parts[2]
    
    # Convert interactor strings into lists of chain IDs
    interactor_1_chain_ids = list(interactor_1)
    interactor_2_chain_ids = list(interactor_2)
    
    # Create output directories if they don't exist
    os.makedirs(args.pdb_out_dir, exist_ok=True)
    os.makedirs(args.rmsd_out_dir, exist_ok=True)
    
    # Always print this message
    print(f"Processing PDB {pdb_id}")
    
    if args.verbose:
        print(f"  Interactor 1: chains {interactor_1_chain_ids}")
        print(f"  Interactor 2: chains {interactor_2_chain_ids}")
    
    # Download the PDB file
    if args.verbose:
        print(f"\nDownloading PDB file...")
    pdbl = PDBList()
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=args.temp_dir, file_format='pdb')
    
    # Load the structure
    if args.verbose:
        print(f"Loading PDB structure from {pdb_file_path}...")
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(pdb_id, pdb_file_path)
    
    # Process interactor 1
    if args.verbose:
        print(f"\nProcessing interactor 1 (chains {interactor_1})...")
    interactor_1_af2_aligned_struct, interactor_1_af2_results = process_af2_alignment(
        structure, pdb_id, interactor_1_chain_ids, verbose=args.verbose
    )
    
    # Process interactor 2
    if args.verbose:
        print(f"\nProcessing interactor 2 (chains {interactor_2})...")
    interactor_2_af2_aligned_struct, interactor_2_af2_results = process_af2_alignment(
        structure, pdb_id, interactor_2_chain_ids, verbose=args.verbose
    )
    
    # Combine RMSD results
    rmsd = {
        'pdb_id': pdb_id,
        'interactor_1': interactor_1_af2_results,
        'interactor_2': interactor_2_af2_results
    }
    
    # Save aligned AF2 structures to PDB files
    # Always print this message
    print(f"Saving aligned structures...")
    interactor_1_af2_aligned_struct_pdb_path = os.path.join(args.pdb_out_dir, f'{pdb_id}_{interactor_1}.pdb')
    interactor_2_af2_aligned_struct_pdb_path = os.path.join(args.pdb_out_dir, f'{pdb_id}_{interactor_2}.pdb')
    
    io1 = PDBIO()
    io1.set_structure(interactor_1_af2_aligned_struct)
    io1.save(interactor_1_af2_aligned_struct_pdb_path)
    # Always print saved files
    print(f"  Saved {interactor_1_af2_aligned_struct_pdb_path}")
    
    io2 = PDBIO()
    io2.set_structure(interactor_2_af2_aligned_struct)
    io2.save(interactor_2_af2_aligned_struct_pdb_path)
    # Always print saved files
    print(f"  Saved {interactor_2_af2_aligned_struct_pdb_path}")
    
    # Save RMSD results to JSON
    rmsd_path = os.path.join(args.rmsd_out_dir, f'{pdb_id}_rmsd.json')
    with open(rmsd_path, 'w') as f:
        json.dump(rmsd, f, indent=2)
    # Always print saved JSON
    print(f"  Saved {rmsd_path}")
    
    if args.verbose:
        print(f"\nRMSD Results:")
        print(json.dumps(rmsd, indent=2))
        print(f"\nProcessing complete!")


if __name__ == '__main__':
    main()

