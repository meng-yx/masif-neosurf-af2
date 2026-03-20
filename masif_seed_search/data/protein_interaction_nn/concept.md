Design a permutation-Invatriant Aggregation architecture on a patch-protein level:
- Use the output of masif-ppi-search as input: 80D descriptors:
    - For each surface vertex, there are 80 encoded values describing the geometric and chemical features of an area within 12A in radius. 
- Goal: design a permutation-Invatriant NN to predict protein-protein interaction on a protein level:
    - Training data:
        - Positives: structures of known protein complexes (p1 chain and p2 chain):
            - p1: use N interface vertices, where N is all the verticies within the interface area. It is expected to vary from protein to protein. 
                - p1_iface_descriptor has a shape of (80, N) (use _straight.npy)
                - Additionally, the real-space coordinate of each vertex is provided 
            - p2: use all (M) vertices, where M is variable and can be either smaller or bigger than N. 
                - p2_descriptor has a shape of (80, M) (use _flipped.npy)
                - Additionally, the real-space coordinate of each vertex is provided 
        - Negatives: randomly sampled PDB structures from unrelated PDB IDs
            - p1: randomly sample N verticies
            - p2: randomly choose a chain from a different PDB structure
    - Idea:
        - The assumption is that if p2 interacts with p1's target patches, there should will be many vertices within p2 with similar descriptors to p1's target patches, and they are oriented in the 3D space such that after applying a transformation to p2, vertices with similar descriptors between p1 and p2 are close to each other. 

