from typing import List
import os
import pandas as pd

def read_full_list(path: str = "full_list.txt"):
    """Read `path` and return a list of non-empty, stripped lines.

    Args:
        path: Path to the file. Defaults to `full_list.txt` in CWD.

    Returns:
        List of strings, one per non-empty line in the file.
    """
    entry_id = []
    asym_id_1 = []
    asym_id_2 = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            entry_id_e, asym_id_1_e, asym_id_2_e = line.strip().split('_')
            entry_id.append(entry_id_e.lower())
            asym_id_1.append(asym_id_1_e)
            asym_id_2.append(asym_id_2_e)
    entries = pd.DataFrame({
        'entry_id': entry_id, 
        'asym_id_1': asym_id_1, 
        'asym_id_2': asym_id_2
    })
    
    return entries
