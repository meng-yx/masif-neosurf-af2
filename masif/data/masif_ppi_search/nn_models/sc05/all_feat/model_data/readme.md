## masif-ppi-search NN model

Retrained on 24.3.2026 - Yanxiang Meng
- Used the original masif repo's training data set (full list of 5902 PDB entries).
- Used updated mixed sampling strategy:
    - Sampling of positive pairs remain unchanged
    - Sampling of negative pairs was updated:
        - Patches from all structures that are not involved in interfaces were pooled, catalogued in cache job, 
        then pairs of patches that do not come from the same PDB accession were randomly sampled to create a 
        hard negative set.
        - Sample 5 times as many pairs of negatives as positives, with weighting of negatives adjusted accordingly
        during training. 
- Saved best model after 102000 iterations. 

