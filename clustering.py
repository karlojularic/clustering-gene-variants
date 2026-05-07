def greedy_clustering(msa, matrix):

    clusters = {}
    for i in range(len(msa)):
        added = False
        for cluster in clusters.values():
            avg_distance = sum(matrix[i][j] for j in cluster) / len(cluster)
            if avg_distance < 0.05:  # Threshold for clustering
                cluster.append(i)
                added = True
                break
        if not added:
            clusters[len(clusters)] = [i]

    return clusters

def cluster_representative(cluster, matrix):
    """
    The representative sequence of a cluster is the medoid: the sequence with the smallest
    average distance to all other sequences in the same cluster
    """
    if len(cluster) == 1:
        return cluster[0]

    best_fragment = None
    best_average_distance = float("inf")
    for fragment_index in cluster:
        total_distance = 0
        for other_fragment_index in cluster:
            if fragment_index != other_fragment_index:
                total_distance += matrix[fragment_index][other_fragment_index]

        average_distance = total_distance / (len(cluster) - 1)
        
        if average_distance < best_average_distance:
            best_average_distance = average_distance
            best_fragment = fragment_index

    return best_fragment