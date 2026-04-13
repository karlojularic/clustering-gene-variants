import numpy as np

def distance_matrix(msa):

    matrix = np.zeros((len(msa), len(msa)), dtype=float)
    print("Calculating distance matrix from MSA...")

    for i in range(len(msa)):
        for j in range(i + 1, len(msa)):
            total = 0
            for a, b in zip(msa[i], msa[j]):
                if a != b:
                    total += 1
            distance = total / len(msa[i])
            matrix[i][j] = distance
            matrix[j][i] = distance  # Symmetric matrix

    return matrix
