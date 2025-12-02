#!/usr/bin/env python3
import subprocess

list_file = "masif/data/masif_ppi_search/lists/small_list.txt"

with open(list_file, "r") as f:
    for line in f:
        entry = line.strip()
        if not entry:
            continue

        print(f"\n=== Running PDB_to_AF2.py for {entry} ===\n")

        subprocess.run(["python", "PDB_to_AF2.py", entry])
