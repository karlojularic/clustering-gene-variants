import numpy as np
import edlib

def distance_matrix(fragments_dict):

    matrix = np.zeros((len(fragments_dict), len(fragments_dict)), dtype=float)
    
    print("Aligning fragments of length 296...")

    for i, frag1 in enumerate(fragments_dict.values()):
        for j, frag2 in enumerate(fragments_dict.values()):
            if i >= j:
                continue

            result = edlib.align(frag1, frag2, task='distance')
            matrix[i][j] = result['editDistance'] / 296
            matrix[j][i] = result['editDistance'] / 296  # Normalize by length of fragments

    return matrix
